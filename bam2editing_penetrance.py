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

def bam2editing_penetrance(out, handle, bam, fasta, minSites, mapQ, baseQ, verbose):
    """Report editing frequency from individual reads"""
    

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
    parser.add_argument("-b", "--bam", nargs=1, 
                        help="input RNA-Seq BAM file")
    parser.add_argument("-f", "--fasta", type=file, 
                        help="reference FASTA file")
    parser.add_argument("-m", "--minSites", default=3,  type=int,
                        help="minimal editing overlapping with read [%(default)s]")
    parser.add_argument("-q", "--mapQ",   default=15, 
                        help="min mapping quality [%(default)s]")
    parser.add_argument("-Q", "--baseQ",  default=20, 
                        help="min base quality [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    bam2editing_penetrance(o.output, o.input, o.bam, o.fasta, o.minSites, o.mapQ, o.baseQ, o.verbose)
    
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
