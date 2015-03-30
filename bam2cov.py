#!/usr/bin/env python
desc="""Report coverage from BAM file. 

By default, the program dumps reads of interest and depth of coverage information. This
speeds up recalculation by the factor of 20X and should take <1% of BAM file size.

Dependencies:
- pysam (sudo easy_install -U pysam)
- bcbio-gff (sudo easy_install -U bcbio-gff)
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 30/03/2015
"""

import os, sys
import pysam, resource
from datetime import datetime
import numpy as np

def parse_bam(bam, dbname, verbose):
    """Return chr2cov"""
    # connect to sqlite3 db
    cnx = sqlite3.connect(dbname)
    cur = cnx.cursor()
    # prepare tables
    for ref in sam.references:
        cur.execute("CREATE TABLE `%s` (start INT, end INT, id INT)"%ref)
        
    #chr2cov[alg.rname][alg.pos]

    cnx.commit()

def load_intervals(fn, verbose):
    """Return chr2intervals and number of entries"""
    chr2intervals = {}
    for i, rec in enumerate(open(fn)):
        # GTF / GFF
        if fn.endswith(('gtf','gff')):
            chrom, source, ftype, s, e, score, strand = rec.split('\t')[:7]
            s, e = int(s)-1, int(e)
        # BED
        else:
            chrom, s, e, name, score, strand = rec.split('\t')[:6]
            s, e = int(s), int(e)
        if strand=="+":
            strand = 1
        else:
            strand = 0
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
        chr2intervals[chrom]=np.array(data, dtype=dtype)        
    return chr2intervals, i
    
def add_reads(blocks, ac, ivals, counts):
    """Add reads"""
    '''# get chromosome name
    c = sam.references[a.rname]
    # count read alignment only once for each block/interval
    selected = c2i[c][np.any([np.any([np.all([s<c2i[c]['start'], c2i[c]['start']<e], axis=0),
                                      np.all([s<c2i[c]['end'], c2i[c]['end']<e], axis=0)], axis=0)
                              for s, e in a.blocks], axis=0)]
    
    d = []
    d+= [np.all([s<ivals['start'], ivals['start']<e], axis=0) for s, e in blocks]
    d+= [np.all([s<ivals['end'], ivals['end']<e], axis=0) for s, e in blocks]
    selected = ivals[np.any(d, axis=0)]'''
    selected = ivals[np.any([np.any([np.all([s<ivals['start'], ivals['start']<e], axis=0),
                                     np.all([s<ivals['end'], ivals['end']<e], axis=0)], axis=0)
                             for s, e in blocks], axis=0)]

    for s, e, strand, entry_id in selected:
        counts[0][entry_id] += 1
    return counts

def _filter(a, mapq=0):
    """Return True if poor quality alignment"""
    if a.mapq<mapq or a.is_secondary or a.is_duplicate or a.is_qcfail:
        return True
            
def buffer_intervals(c2i, ivals, sam, a, maxp, pref, bufferSize):
    """Return invervals buffer for faster selection"""
    if a.aend>maxp or a.rname != pref:
        # get ref chrom
        c = sam.references[a.rname]
        s, e = maxp, maxp+bufferSize
        # update intervals
        #sys.stderr.write(" buffering %s:%s-%s ...\r"%(c,s,e))
        if c in c2i:
            ivals = c2i[c][np.any([np.all([s<c2i[c]['start'], c2i[c]['start']<e], axis=0),
                                   np.all([s<c2i[c]['end'], c2i[c]['end']<e], axis=0)], axis=0)]
        else:
            ivals = np.empty_like(ivals)
        #sys.stderr.write("  %s entries loaded!\r"%len(ivals))
        # store current reference and max position
        pref = a.rname
        maxp = e
    return ivals, maxp, pref

def parse_bam(bam, mapq, c2i, entries, bufferSize, verbose):
    """Parse BAM and return counts for sense/antisense of each interval"""
    counts = (np.zeros(entries, dtype='uint16'), np.zeros(entries, dtype='uint16'))
    # open BAM
    sam = pysam.AlignmentFile(bam)
    # keep info about previous read 
    pa, ac = 0, 0
    # keep info about intervals, max position and current reference
    ivals, maxp, pref = [], 0, 0
    for i, a in enumerate(sam, 1):
        if not i%1e5:
            sys.stderr.write(' %i\r'%(i, )) #resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
        #if i>1e6: break
        # filter poor quality
        if _filter(a, mapq):
            continue
        if not pa:
            pa, ac = a, 1
            continue
        # check if similar to previous
        if pa.pos==a.pos and pa.cigarstring==a.cigarstring:
            ac += 1
        else:
            ivals, maxp, pref = buffer_intervals(c2i, ivals, sam, pa, maxp, pref, bufferSize)
            counts = add_reads(pa.blocks, ac, ivals, counts)
            pa, ac = a, 1
    # add last alignment
    if not _filter(a, mapq):
        counts = add_reads(a.blocks, ac, ivals, counts)
    return counts
    
def bam2cov(bam, bed, out=sys.stdout, mapq=0, bufferSize=1000000, verbose=1):
    """Calculate coverage for genome intervals."""
    # load intervals
    if verbose:
        sys.stderr.write("Loading intervals...\n")
    c2i, entries = load_intervals(bed, verbose)
    # parse alignments
    if verbose:
        sys.stderr.write("Parsing alignments...\n")
    counts = parse_bam(bam, mapq, c2i, entries, bufferSize, verbose)
        
    print counts[0].sum()
    
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
                        help="BED/GTF/GFF file")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-p", "--ploidy",    default=2, type=int, 
                        help="ploidy          [%(default)s]")
    parser.add_argument("-q", "--mapq",      default=20, type=int, 
                        help="min mapping quality for variants [%(default)s]")
    parser.add_argument("--bufferSize",      default=1e6,  type=int, 
                        help="buffer size for intervals [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    #initialise structural variants
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
