#!/usr/bin/env python
desc="""Report coverage from BAM file. 
Support spliced alignments and mapping quality filtering.
By default ignores secondary alignments, duplicates and quality failed reads. 

bam2cov_pool.py is nearly 2x faster on large BAM files than bam2cov.py. 

TDB:
- define minimum overlap

CHANGELOG:
v1.1
- skip counting read2 for paired-end libs (like samtools view -F 128)
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 30/03/2015
"""

import os, sys, pysam, re, subprocess
from datetime import datetime
import numpy as np
from multiprocessing import Pool

cigarPat = re.compile('\d+\w')

def load_intervals(fn, verbose):
    """Return chr2intervals and number of entries"""
    chr2intervals = {}
    for i, rec in enumerate(open(fn)):
        if rec.startswith('#') or not rec.strip():
            continue
        rdata = rec.split('\t')
        score, strand = 1, "+"    
        # GTF / GFF
        if fn.endswith(('gtf','gff')):
            chrom, source, ftype, s, e, score, strand = rdata[:7]
            s, e = int(s)-1, int(e)
        # BED
        else:
            # unstranded intervals
            if len(rdata)<6:
                chrom, s, e = rdata[:3]
            else:
                chrom, s, e, name, score, strand = rdata[:6]
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
    """Return True if poor quality alignment or second in pair."""
    if a.mapq<mapq or a.is_secondary or a.is_duplicate or a.is_qcfail or alg.is_read2:
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
    d  = [np.all([ (s+e)/2.>=ivals['start'], (s+e)/2.<=ivals['end'] ], axis=0) for s, e in blocks]
    #d  = [np.all([ s>=ivals['start'], s<=ivals['end'] ], axis=0) for s, e in blocks]
    # end overlapping with interval
    #d += [np.all([ e>=ivals['start'], e<=ivals['end'] ], axis=0) for s, e in blocks]
    # interval inside read
    #d += [np.all([ s< ivals['start'], e> ivals['end'] ], axis=0) for s, e in blocks]
    # select intervals fulfilling any of above
    selected = ivals[np.any(d, axis=0)]
    cminus = strands.count(True)
    cplus  = len(strands)-cminus #strands.count(False)   
    return cminus, cplus, list(selected)

def cigar2blocks(pos, cigar):
    """Return alignment blocks and aend. cigar2block doesn't break
    block upon insertion (I), while pysam does."""
    blocks = [[pos, pos]]
    for c in cigarPat.findall(cigar):
        bases, s = int(c[:-1]), c[-1]
        '''if s in ('S','H','I'):
            continue
        # match'''
        if s == 'M':
            blocks[-1][-1] += bases
        # intron
        elif s in ("N","D"):
            pos = blocks[-1][-1] + bases
            blocks.append([pos, pos])
        '''# insertion - skipped earlier
        elif s =="I":
            pos = blocks[-1][-1] # + bases
            blocks.append([pos, pos])'''
    #for i in range(len(blocks)): blocks[i] = tuple(blocks[i])
    return blocks

def is_reverse(flag):
    """Return true if read aligned to reverse strand"""
    if int(flag) & 16:
        return True
    return False
    
def alignment_iterator_samtools(bam, mapq, verbose):
    """Iterate aligments from BAM using samtools view subprocess"""
    # start samtools view subprocess
    cmd0  = ["samtools", "view", "-q%s"%mapq, "-F3968", bam]
    # 3968 = skip: read2, secondary, QC fail, duplicates and supplementary (http://broadinstitute.github.io/picard/explain-flags.html)
    proc0 = subprocess.Popen(cmd0, bufsize=-1, stdout=subprocess.PIPE)
    out0  = proc0.stdout
    # start cut subprocess
    cmd1  = ["cut", "-f2-4,6"]
    proc1 = subprocess.Popen(cmd1, bufsize=-1, stdout=subprocess.PIPE, stdin=out0)
    out1  = proc1.stdout
    if verbose:
        sys.stderr.write(" "+" | ".join([" ".join(cmd0), " ".join(cmd1)+"\n"]))
    # keep info about previous read 
    strands = []
    # iterate only ok alignments
    for i, line in enumerate(out1, 1):
        #if i>1e5: break
        if verbose and not i%1e5:
            sys.stderr.write(' %i \r'%i)
        flag, rname, pos, cigar = line[:-1].split('\t')
        pos = int(pos)
        # get blocks, aend & strand
        blocks  = cigar2blocks(pos, cigar)
        aend    = blocks[-1][-1]
        reverse = is_reverse(flag)
        if not strands:
            ppos, pcigar, paend, pblocks, prname = pos, cigar, aend, blocks, rname
            strands = [reverse]
            continue
        # check if similar to previous
        if ppos==pos and pcigar==cigar:
            strands.append(reverse)
        else:
            yield pblocks, strands, ppos, paend, prname
            # store current entry
            ppos, pcigar, paend, pblocks, prname = pos, cigar, aend, blocks, rname
            strands = [reverse]
    # info
    if verbose:
        sys.stderr.write(' %i alignments processed.\n'%i)
    yield pblocks, strands, ppos, paend, prname
        
def parse_bam(bam, mapq, c2i, entries, bufferSize, verbose, nprocs=4):
    """Parse BAM and return counts for sense/antisense of each interval"""
    counts = ([0]*entries, [0]*entries)
    # init workers
    p = Pool(nprocs, initializer=init_args, initargs=(c2i, bufferSize))
    # init iterator
    iterator = alignment_iterator_samtools(bam, mapq, verbose) 
    # parse results
    for cminus, cplus, selected in p.imap_unordered(worker, iterator, chunksize=10000):
        # store info
        for s, e, reverse, ival in selected:
            # - transcript on reverse
            if reverse:
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
