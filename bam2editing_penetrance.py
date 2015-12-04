#!/usr/bin/env python
desc="""Report frequency of editing from individual reads. 
 
CHANGELOG:
v1.0
- 

TBD:
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw, 2/12/2015
"""

import argparse, os, re, sys
from datetime import datetime
import subprocess
import numpy as np
import pysam

def bed2intervals(handle, window=1000):
    """Yield collapsed RNA editing cases in given window"""
    data = []
    for l in handle:
        chrom, s, e, change, score, strand, freq = l[:-1].split('\t')[:7]
        s, e = map(int, (s, e))
        # store only if same chrom and distance between prev e and curr start < window
        if not data:
            pass
        elif data[-1][0]!=chrom or data[-1][2]-s > window:
            yield data
            data = []
        data.append((chrom, s, e, change))
    yield data

def bam2editing_penetrance(out, handle, bam, minSites, mapQ, baseQ, verbose):
    """Report editing frequency from individual reads"""
    # NOTE USING SMALL WINDOW WILL LIKELY SKIP EXON-INTRON SPANNING READS! OR THE WILL BE COUNTED TWICE
    sam = pysam.Samfile(bam)
    for intervals in bed2intervals(handle):
        if len(intervals)<minSites:
            continue
        # get edited ref positions A>G to ['A', 'G']
        pos2info = { s: change.split('>') for chrom, s, e, change in intervals } 
        # get region
        chrom, s = intervals[0][:2]
        e = intervals[-1][2]
        # fetch reads and parse
        for r in sam.fetch(chrom, s, e):
            # skip low quality reads
            if r.mapq<mapQ or r.is_secondary or r.is_duplicate or r.is_qcfail:
                continue
            # process aligned pairs of bases
            i = k = 0
            for qi, ri in r.get_aligned_pairs():
                # skip q dels or if not edited
                if not qi or ri not in pos2info:
                    continue
                # skip low base call quality
                if r.query_qualities[qi]<baseQ:
                    continue
                i += 1
                ref, alt = pos2info[ri]
                if   r.seq[qi] == ref:
                    i += 1
                elif r.seq[qi] == alt:
                    i += 1
                    k += 1
            if i<minSites:
                continue
            out.write("%s:%s-%s\t%s\t%s\t%s\t%s\n"%(chrom, s, e, r.qname, k, i, 1.0*k/i))

def main():

    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0a')
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-i", "--input", type=file, default=sys.stdin,  
                        help="input stream with RNA editing positions [stdin]")
    parser.add_argument("-b", "--bam", required=True, 
                        help="input RNA-Seq BAM file")
    #parser.add_argument("-f", "--fasta", type=file, 
    #                    help="reference FASTA file")
    parser.add_argument("-m", "--minSites", default=3,  type=int,
                        help="minimal editing overlapping with read [%(default)s]")
    parser.add_argument("-q", "--mapQ",   default=15, 
                        help="min mapping quality [%(default)s]")
    parser.add_argument("-Q", "--baseQ",  default=20, 
                        help="min base quality [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    bam2editing_penetrance(o.output, o.input, o.bam, o.minSites, o.mapQ, o.baseQ, o.verbose)
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
