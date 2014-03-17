#!/usr/bin/env python
desc="""Generate batch file for IGV
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow, 17/03/2014
"""

import os, sys
import pickle, pysam, resource
from datetime import datetime
import numpy as np
from scipy import stats

def bed2batch(bed, out, session, outdir, ext, offset, verbose):
    """Generates IGV batch script."""
    if session:
        out.write("new\nload %s\n"%session)
    if outdir:
        out.write("snapshotDirectory %s\n"%outdir)

    #process bed
    for line in bed:
        if line.startswith("#"):
            continue
        chrom, s, e, name, score = line.split('\t')[:5]
        s, e = int(s), int(e)
        soff, eoff = s-500, e+500
        if soff<1:
            soff = 1
        coords = "%s:%s-%s" % (chrom, soff, eoff)            
        outfn = "%s.%s:%s-%s.%s" % (name, chrom, s, e, ext)
        out.write("goto %s\nsnapshot %s\n" % (coords, outfn))

        
def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-i", "--bed",       default=sys.stdin, type=file, 
                        help="BED file        [stdin]")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-s", "--session",   default="", required=True, 
                        help="session file")
    parser.add_argument("--outdir",          default="", required=True, 
                        help="output dir")
    parser.add_argument("--ext",             default="png",
                        choices=['png', 'jpg', 'svg'],  
                        help="snapshots format [%(default)s]")
    parser.add_argument("--offset",          default=500, type=int, 
                        help="start/end offset [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    bed2batch(o.bed, o.output, o.session, o.outdir, o.ext, o.offset, o.verbose)

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
