#!/usr/bin/env python
desc="""Identify enhanvers (eRNAs) from RNA-Seq.
Enhancers are defined as genomic regions having similar level of expression. 

CHANGELOG:
v1.1
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw, 8/10/2015
"""

import argparse, os, pysam, sys
from datetime import datetime
import numpy as np

def parse_algs(bamFn, mapq, verbose):
    """Parse alignements."""
    cov = rname = None
    # parse algs
    bam = pysam.Samfile(bamFn)
    for alg in bam: #bam.fetch(reference='1'):
        # skip low mapq, unmapped and secondary
        if alg.mapq < mapq or alg.rname<0 or alg.is_secondary:
            continue
        if alg.rname != rname:
            # report previous
            if cov:
                yield bam.references[rname], cov
            # reset rname and chr2cov
            rname = alg.rname
            # stores +/- read counts
            cov = [np.zeros(bam.lengths[rname]), np.zeros(bam.lengths[rname])]
        # get strand
        strand = 0
        if alg.is_reverse:
            strand = 1
        # update coverage
        for s, e in alg.blocks:
            cov[strand][s:e] += 1
    # report last
    if cov:
        yield bam.references[rname], cov

def cov2balanced(cov, minReads, delta, verbose):
    """Return intervals with small delta"""
    s = e = previous = 0
    for i, (a, b) in enumerate(zip(cov[0], cov[1])):
        # break interval
        if a<minReads or b<minReads or abs(np.log(1.0*a/b)) > delta:
            if previous:
                e = i
                yield s, e, cov[0][s:e], cov[1][s:e]
            previous = 0
        # start new interval
        elif not previous:
            previous = 1
            s = i
            
    if previous:
        e = i
        yield s, e, cov[0][s:e], cov[1][s:e]
        
def bam2enhancers(out, bam, mapq, minLength, minReads, delta, verbose):
    """Report enhancer sequences"""
    for i, (chrom, cov) in enumerate(parse_algs(bam, mapq, verbose), 1):
        for s, e, plus, minus in cov2balanced(cov, minReads, delta, verbose):
            if e-s<minLength:
                continue
            log2 = np.log2(1.0*sum(plus)/sum(minus))
            out.write("%s\t%s\t%s\t%i\t%i\t%.3f\n"%(chrom, s, e, np.median(plus), np.median(minus), log2))
        
def main():

    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.11')
    parser.add_argument("-i", "--bam", type=file, 
                        help="input BAM file")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("w"), 
                        help="output stream [stdout]")
    parser.add_argument("-q", "--mapq",      default=10, type=int, 
                        help="min mapping quality for variants [%(default)s]")
    parser.add_argument("-l", "--minLength",  default=25, type=int, 
                        help="min number of reads [%(default)s]")
    parser.add_argument("-m", "--minReads",  default=10, type=int, 
                        help="min number of reads [%(default)s]")
    parser.add_argument("-d", "--delta",     default=0.1, type=int, 
                        help="max difference between +/- strand [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
    
    bam2enhancers(o.output, o.bam.name, o.mapq, o.minLength, o.minReads, o.delta, o.verbose)
    
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
