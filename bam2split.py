#!/usr/bin/env python3
desc="""Split primary alignments from BAM file into N files.
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 16/07/2021
"""

import os, sys, pysam
from datetime import datetime

def is_qcfail(a, mapq=20):
    """Return True if read fails QC"""
    # we don't mind supplementary algs in general
    if a.mapq<mapq or a.is_qcfail or a.is_secondary or a.is_duplicate: return True

def bam2split(outbase, bam, nfiles, mapq):
    """Split primary alignments passing mapping quality threshold into N output files"""
    outdir = os.path.dirname(outbase)
    if not os.path.isdir(outdir): os.makedirs(outdir)
    # get input sam
    sam = pysam.AlignmentFile(bam)
    # open output files
    out = [pysam.AlignmentFile("%s.%s.bam"%(outbase, n), "wb", header=sam.header)
           for n in range(1, nfiles+1)]
    # split alignments
    i = 0
    for ai, a in enumerate(sam, 1):
        if not ai%1000: sys.stderr.write(" %s / %s   \r"%(ai, sam.mapped))
        if is_qcfail(a, mapq): continue
        out[i%nfiles].write(a)
        i += 1
    sys.stderr.write("%s out of %s alignments written to %s output files\n"%(i, sam.mapped, nfiles))
    # close & index
    for o in out:
        o.close()
        pysam.index(o.filename)

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0a')
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-i", "--bam", required=True, help="input BAM file")
    parser.add_argument("-n", "--nfiles", default=2, type=int, help="number of files to split into [%(default)s]")
    parser.add_argument("-o", "--outbase", required=True, help="output base name")
    parser.add_argument("--mapq", default=15, type=int, help="min mapping quality [%(default)s]")

    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    bam2split(o.outbase, o.bam, o.nfiles, o.mapq)
        
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
