#!/usr/bin/env python
desc="""Scan alignments (SAM/BAM) for structural variants.

Detection of deletions, duplications and inversions from paired reads is implemented.
In addition, deletions and duplications are detected from deviations from mean depth
of coverage.

By default, the program dumps reads of interest and depth of coverage information. This
speeds up recalculation by the factor of 20X and should take <1% of BAM file size.

To be done:
+ SV detection
-- insertions testing
-- translocations
-- inversions
+ split read alignment
+ rlen learning
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 13/03/2014
"""

import os, sys
from datetime import datetime
import numpy as np
from scipy import stats, signal

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
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-p", "--ploidy",    default=2, type=int, 
                        help="ploidy          [%(default)s]")
    parser.add_argument("-q", "--mapq",      default=20, type=int, 
                        help="min mapping quality for variants [%(default)s]")
    parser.add_argument("--rlen",            default=None, type=int, 
                        help="read length     [get from data]")
    parser.add_argument("-c", "--covD",      default=0.33, type=float, 
                        help="min coverage change to call deletion/duplication [%(default)s]")
    parser.add_argument("--cov_frac",        default=0.1, type=float, 
                        help="min fraction of local depth to call variation [%(default)s]")
    parser.add_argument("--dup_isize_frac",  default=0.9, type=float, 
                        help="min duplication size as insert size fraction [%(default)s]")
    parser.add_argument("--cnv_size",        default=1000, type=int, 
                        help="min CNV size from depth of coverage [%(default)s]")
    parser.add_argument("--merge",           default=False,  action="store_true",
                        help="merge read pairs variants using depth of coverage variants")
    parser.add_argument("--nodump",          default=False,  action="store_true",
                        help="dump SV reads for faster recalculations")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    #initialise structural variants
    sv = SVs(o.bam, out=o.output, mapq=o.mapq, ploidy=o.ploidy, covD=o.covD, \
             cov_frac=o.cov_frac, rlen=o.rlen, dup_isize_frac=o.dup_isize_frac, \
             cnv_size=o.cnv_size, merge=o.merge, \
             nodump=o.nodump, verbose=o.verbose)
    #call variants in all chromosomes
    sv.parse()

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
