#!/usr/bin/env python
desc="""Calculate completeness of transcripts

Assumes sense read orientation (as in direct RNA-seq ONT). 
Support spliced alignments and mapping quality filtering.
By default ignores secondary alignments, duplicates and quality failed reads.
Takes into account only the first read from pair. 

USAGE: f=minimap2.ref10/20190219.bam; ../src/bam2transcript_completeness.py -v -i $f -b ../ref10/DANRE.gff3 | plot_hist.py -b 25 -t $f -o $f.transcript_completeness.png

"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Oxford, 27/02/2019
"""

import os, sys, pysam
from datetime import datetime
import numpy as np

def description2name(desc):
    """Return name from GFF/GTF description"""
    k2v = {kv.split('=')[0]: "=".join(kv.split('=')[1:])
           for kv in desc.split(';') if kv}
    if "ID" in k2v:
        return k2v["ID"]
    else:
        return "NoName"

def load_intervals(fn, verbose=1, ftypefilter="gene"):
    """Return chr2intervals and number of entries"""
    chr2intervals = {}
    for i, rec in enumerate(open(fn)):
        if rec.startswith('#') or not rec.strip():
            continue
        rdata = rec.split('\t')
        score, strand = 1, "+"    
        # GTF / GFF
        if fn.endswith(('gtf', 'gff', 'gff3')):
            chrom, source, ftype, s, e, score, strand = rdata[:7]
            s, e = int(s)-1, int(e)
            name = description2name(rdata[8])
            if name.startswith('gene:'): name=name[5:]
        # BED
        else:
            # unstranded intervals
            if len(rdata)<6:
                chrom, s, e = rdata[:3]
                name = i
            else:
                chrom, s, e, name, score, strand = rdata[:6]
            s, e = int(s), int(e)
            ftype = ''
        # filter by feature type
        if ftypefilter and ftypefilter!=ftype:
            continue
        if strand=="+":
            strand = 0
        else:
            strand = 1
        # add chromosome
        if chrom not in chr2intervals:
            chr2intervals[chrom] = []
        # store interval
        data = (s, e, strand, name)
        chr2intervals[chrom].append(data)

    # define numpy datatype
    dtype = np.dtype({'names':   ['start',  'end',    'strand', 'entry_id'], \
                      'formats': ['uint32', 'uint32', 'bool_', 'object']})
    for chrom, data in chr2intervals.iteritems():
        chr2intervals[chrom] = np.array(data, dtype=dtype)
    return chr2intervals, i
    
def _filter(a, mapq=0):
    """Return True if poor quality alignment"""
    if a.mapq<mapq or a.is_secondary or a.is_duplicate or a.is_qcfail:
        return True
            
def buffer_intervals(c2i, ivals, rname, pos, aend, maxp, pref, bufferSize):
    """Return invervals buffer for faster selection"""
    if rname != pref:
        maxp = 0
    if aend > maxp:
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
        pref = rname
        maxp = e
    #print(ivals, rname, pos, aend, pref, bufferSize, c2i.keys)
    return ivals, maxp, pref

def count_overlapping_intervals(blocks, is_reverse, ivals, out, verbose=0):
    """Count overlapping intervals with given read alignment.
    The algorithm support spliced alignments. """
    # skip if not ivals
    if not len(ivals):
        return 0.
    ## get intervals overlapping with given alignment blocks
    # start overlapping with interval - here I'm assuming sense read orientation
    d  = [np.all([ (s+e)/2.>=ivals['start'], (s+e)/2.<=ivals['end'], is_reverse==ivals['strand'] ], axis=0) for s, e in blocks]
    # select intervals fulfilling any of above
    selected = ivals[np.any(d, axis=0)]#; print selected
    # check if any matches, as sometimes empty cause problems
    name, overlap = "-", 0.0
    for s, e, strand, _name in selected:
        # get only parts of the read the aligned to gene/transcript
        re = blocks[-1][-1] if blocks[-1][-1] < e else e
        rs = blocks[0][0] if blocks[0][0] > s else s
        # get best match
        _overlap = 1.*(re-rs)/(e-s)
        if _overlap>overlap:
            name, overlap = _name, _overlap
    #print selected, is_reverse, blocks[0][0], blocks[-1][-1], overlap        
    out.write("%s\t%s\n"%(name, overlap))
            
def alignment_iterator(bam, mapq, verbose):
    """Iterate alignments from BAM"""
    # open BAM
    sam = pysam.AlignmentFile(bam)
    # count alg quality ok
    qok = 0
    # keep info about previous read 
    pa, strands = 0, []
    for i, a in enumerate(sam, 1):
        #if i>1e5: break
        #if i<84*1e5: continue
        if verbose and not i%1e5:
            sys.stderr.write(' %i algs; %i ok \r'%(i, qok))
        # filter poor quality
        if _filter(a, mapq):
            continue
        qok += 1
        yield a.blocks, a.is_reverse, a.pos, a.aend, sam.references[a.rname]
    # info
    if verbose:
        sys.stderr.write(' %i alignments processed.\n'%i)
    
def bam2overlap(bam, bed, out=sys.stdout, mapq=0, bufferSize=1000000, verbose=1):
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
    # write header
    out.write("# gene\tread overlap\n")
    ivals, maxp, pref = [], 0, 0
    for blocks, strands, pos, aend, rname in alignment_iterator(bam, mapq, verbose): #alignment_iterator_samtools(bam, mapq, verbose):
        # update ivals
        ivals, maxp, pref = buffer_intervals(c2i, ivals, rname, pos, aend, maxp, pref, bufferSize)
        # add alignments
        count_overlapping_intervals(blocks, strands, ivals, out, verbose)
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.1')   
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
    bam2overlap(o.bam, o.bed, o.output, o.mapq, o.bufferSize, o.verbose)
 
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
