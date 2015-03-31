#!/usr/bin/env python
desc="""Report coverage from BAM file. 
Support spliced alignments and mapping quality filtering.
By default ignores secondary alignments, duplicates and quality failed reads. 

Dependencies:
- pysam (sudo easy_install -U pysam)

TDB:
- define minimum overlap
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 30/03/2015
"""

import os, sys, pysam
from datetime import datetime
import numpy as np
from multiprocessing import Pool, Process, Queue

def load_intervals(fn, verbose):
    """Return chr2intervals and number of entries"""
    chr2intervals = {}
    for i, rec in enumerate(open(fn)):
        if rec.startswith('#') or not rec.strip():
            continue
        # GTF / GFF
        if fn.endswith(('gtf','gff')):
            chrom, source, ftype, s, e, score, strand = rec.split('\t')[:7]
            s, e = int(s)-1, int(e)
        # BED
        else:
            chrom, s, e, name, score, strand = rec.split('\t')[:6]
            s, e = int(s), int(e)
        if strand=="+":
            strand = 0
        else:
            strand = 1
        # add chromosome
        if chrom not in chr2intervals:
            chr2intervals[chrom] = []
        # store interval
        data = (s, e, strand, i)
        chr2intervals[chrom].append(data)

    # define numpy datatype
    dtype = np.dtype({'names':   ['start',  'end',    'strand', 'entry_id'], \
                      'formats': ['uint32', 'uint32', 'bool_', 'uint32']})
    for chrom, data in chr2intervals.iteritems():
        chr2intervals[chrom] = np.array(data, dtype=dtype)
    return chr2intervals, i+1
    
def _filter(a, mapq=0):
    """Return True if poor quality alignment"""
    if a.mapq<mapq or a.is_secondary or a.is_duplicate or a.is_qcfail:
        return True
            
def buffer_intervals(c2i, ivals, pos, aend, rname, maxp, pref, bufferSize):
    """Return invervals buffer for faster selection"""
    if rname != pref:
        maxp = 0
    if aend>maxp:
        # get ref chrom
        c = rname
        s, e = pos, aend+bufferSize
        # update intervals
        if c in c2i:
            # select intervals that either start, end or encompass current window/buffer
            ivals = c2i[c][np.any([np.all([ c2i[c]['start']>=s, c2i[c]['start']<=e ], axis=0),
                                   np.all([ c2i[c]['end']  >=s, c2i[c]['end']  <=e ], axis=0),
                                   np.all([ c2i[c]['start']< s, c2i[c]['end']  > e ], axis=0)], axis=0)]
        else:
            ivals = [] 
        #sys.stderr.write(" new buffer with %s intervals: %s:%s-%s\n"%(len(ivals),c,s,e))
        # store current reference and max position
        pref = c
        maxp = e
    return ivals, maxp, pref

def init_args(*args):
    global c2i, bufferSize, ivals, maxp, pref
    c2i, bufferSize = args
    # keep info about intervals, max position and current reference
    ivals, maxp, pref = [], 0, 0
    
def worker(args):
    """Count overlapping intervals with given read alignment.
    The algorithm support spliced alignments. """
    global c2i, bufferSize, ivals, maxp, pref
    # read args
    blocks, strands, pos, aend, rname = args
    # buffer intervals
    ivals, maxp, pref = buffer_intervals(c2i, ivals, pos, aend, rname, maxp, pref, bufferSize)
    # skip if not ivals
    if not len(ivals):
        return 0, 0, []
    ## get intervals overlapping with given alignment blocks
    # start overlapping with interval
    d  = [np.all([ s>=ivals['start'], s<=ivals['end'] ], axis=0) for s, e in blocks]
    # end overlapping with interval
    d += [np.all([ e>=ivals['start'], e<=ivals['end'] ], axis=0) for s, e in blocks]
    # interval inside read
    d += [np.all([ s< ivals['start'], e> ivals['end'] ], axis=0) for s, e in blocks]
    # select intervals fulfilling any of above
    selected = ivals[np.any(d, axis=0)]
    cminus, cplus = strands.count(True), strands.count(False)    
    return cminus, cplus, list(selected)

def alignment_iterator(inQ, bam, mapq, nprocs, verbose):
    """Yield merged alignments that overlap entirely.
    This means start position and cigar is the same."""
    # open BAM
    sam = pysam.AlignmentFile(bam)
    # count alg quality ok
    qok = 0
    # keep info about previous read 
    pa, strands = 0, []
    for i, a in enumerate(sam, 1):
        #if i>3*1e5: break
        #if i<84*1e5: continue
        if verbose and not i%1e5:
            sys.stderr.write(' %i algs; %i ok \r'%(i, qok))
        # filter poor quality
        if _filter(a, mapq):
            continue
        qok += 1
        if not pa:
            pa = a
            continue
        # check if similar to previous
        if pa.pos==a.pos and pa.cigarstring==a.cigarstring:
            strands.append(a.is_reverse)
        else:
            #yield pa.blocks, strands, pa.pos, pa.aend, sam.references[pa.rname]
            inQ.put((pa.blocks, strands, pa.pos, pa.aend, sam.references[pa.rname]))
            # store current entry
            pa = a
            strands = [a.is_reverse]
    if verbose:
        sys.stderr.write(' %i alignments processed.\n'%i)
    inQ.put((pa.blocks, strands, pa.pos, pa.aend, sam.references[pa.rname]))
    for p in range(nprocs):
        inQ.put('STOP')
    
def parse_bam(bam, mapq, c2i, entries, bufferSize, verbose, nprocs=4):
    """Parse BAM and return counts for sense/antisense of each interval"""
    counts = (np.zeros(entries, dtype='uint32'), np.zeros(entries, dtype='uint32'))
    #get in & out queue
    inQ = Queue()
    # start job populating input queue           
    Process(target=alignment_iterator, args=(inQ, bam, mapq, nprocs, verbose)).start()
    # init workers
    p = Pool(nprocs, initializer=init_args, initargs=(c2i, bufferSize))
    # init iterator
    iterator = iter(inQ.get, 'STOP') #alignment_iterator(bam, mapq, verbose)
    # parse results
    for cminus, cplus, selected in p.imap_unordered(worker, iterator, chunksize=10000):
        # store info
        for s, e, strand, ival in selected:
            # - transcript on reverse
            if strand:
                counts[0][ival] += cminus
                counts[1][ival] += cplus
                # + transcript on forward
            else:
                counts[0][ival] += cplus
                counts[1][ival] += cminus    
    return counts
    
def bam2cov(bam, bed, out=sys.stdout, mapq=0, bufferSize=1000000, verbose=1):
    """Calculate coverage for genome intervals."""
    # load intervals
    if verbose:
        sys.stderr.write("Loading intervals...\n")
    c2i, entries = load_intervals(bed, verbose)
    if verbose:
        sys.stderr.write(" %s intervals from %s chromosomes loaded!\n"%(entries, len(c2i)) )
    # parse alignments & count interval overlaps
    if verbose:
        sys.stderr.write("Parsing alignments...\n")
    counts = parse_bam(bam, mapq, c2i, entries, bufferSize, verbose)
    if verbose:
        sys.stderr.write(" sense / antisense alignments: %s / %s\n" % (sum(counts[0]), sum(counts[1])) )
    # report
    #return 
    for sense, antisense, line in zip(counts[0], counts[1], open(bed)):
        out.write("%s\t%s\t%s\n"%(line[:-1], sense, antisense))
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--bam", required=True,       
                        help="BAM file")
    parser.add_argument("-b", "--bed", required=True,       
                        help="BED/GTF/GFF interval file")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-q", "--mapq",      default=10, type=int, 
                        help="min mapping quality for variants [%(default)s]")
    parser.add_argument("--bufferSize",      default=100000,  type=int, 
                        help="buffer size for intervals [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    # calculate coverage from bam for given intervals
    bam2cov(o.bam, o.bed, o.output, o.mapq, o.bufferSize, o.verbose)
 
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
