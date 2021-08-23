#!/usr/bin/env python2
desc="""Generate batch file for IGV
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Mizerow, 17/03/2014
"""

import os, sys
import pickle, pysam, resource
from datetime import datetime
import numpy as np
from scipy import stats

def bed2batch(bed, out, session, bam, genome, outdir, ext, offset, noname, \
              replace, verbose):
    """Generates IGV batch script."""
    init = "new\n"
    if session:
        init += "load %s\n"%session
    if bam:
        init += "load %s\n"%",".join(bam)
    if genome:
        init += "genome %s\n"%genome
    if outdir:
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        init += "snapshotDirectory %s\n"%outdir
    init += "maxPanelHeight 150\ncollapse\nexpand Gene\n"
    out.write(init)
    #process bed
    i = 0
    for line in bed:
        if line.startswith("#"):
            continue
        chrom, s, e, name = line[:-1].split('\t')[:4]
        s, e = int(s), int(e)
        size = e-s
        #show 2*size window
        soff, eoff = int(s-0.5*size), int(e+0.5*size)
        #increase window size
        if eoff - soff < 2*offset:
            soff, eoff = s-offset, e+offset
        if soff<1:
            soff = 1
        coords = "%s:%s-%s" % (chrom, soff, eoff)
        if noname:
            outfn = "%s:%s-%s.%s" % (chrom, s, e, ext)            
        else:
            outfn = "%s:%s-%s.%s.%s" % (chrom, s, e, name[:100], ext)
        outpath = os.path.join(outdir, outfn)
        if not replace and os.path.isfile(outpath):
            continue
        outline = "goto %s\nsnapshot %s\n" % (coords, outfn)
        out.write(outline)
        #write 3x first event, so size is adjusted
        if not i:
            out.write(outline+outline)
        i += 1
    out.write("exit")
        
def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-v", "--verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-r", "--replace",  default=False, action="store_true", help="replace files")    
    parser.add_argument("-i", "--bed",       default=sys.stdin, type=file, 
                        help="BED file        [stdin]")
    parser.add_argument("--bam",             default=[], nargs="*",
                        help="BAM files")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-s", "--session",   default="", 
                        help="session file")
    parser.add_argument("-g", "--genome",    default="", 
                        help="genome id")
    parser.add_argument("--outdir",          default="", required=True, 
                        help="output dir")
    parser.add_argument("--ext",             default="png",
                        choices=['png', 'jpg', 'svg'],  
                        help="snapshots format [%(default)s]")
    parser.add_argument("--offset",          default=500, type=int, 
                        help="min start/end offset [%(default)s]")
    parser.add_argument('--noname', default=False, action='store_true',
                        help="skip BED name columns in output filename")   
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    if not o.session and not o.bam:
        sys.stderr.write("BAM file or session file need to be provided!\n")
        sys.exit(1)
    bed2batch(o.bed, o.output, o.session, o.bam, o.genome, o.outdir, o.ext, \
              o.offset, o.noname, o.replace, o.verbose)

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
