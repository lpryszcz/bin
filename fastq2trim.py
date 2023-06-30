#!/usr/bin/env python
desc="""Trim fastq alignments.
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Bratislava, 16/09/2015
"""

import os, sys
from datetime import datetime
from Bio import SeqIO

def fastq2trim(handle, out, start, end, verbose):
    """Parse fastQ and report trimmed sequences"""
    for r in SeqIO.parse(handle, "fastq"):
        out.write(r[start:end].format("fastq"))

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--input", default=sys.stdin, type=argparse.FileType('rt'), 
                        help="input stream    [stdin]")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType('wt'), 
                        help="output stream   [stdout]")
    parser.add_argument("-s", "--start",  default=0, type=int, 
                        help="read start      [%(default)s]")
    parser.add_argument("-e", "--end",   default=10000000, type=int, 
                        help="read end        [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    fastq2trim(o.input, o.output, o.start, o.end, o.verbose)
    
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
