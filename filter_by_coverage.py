#!/usr/bin/env python2
desc="""Filter BED entries by coverage
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com
Mizerow, 2014/02/14
"""

import argparse, gzip, os, sys
import pysam
import numpy as np
from datetime import datetime

def filter_by_coverage(bamfn, bedstream, covfr, verbose):
    """Report only windows with at least mincov of mean"""
    #calculate mean cov
    sam = pysam.Samfile(bamfn)

    #mean covevarage  = aligned reads / genome size
    meancov = 1.0 * sam.mapped / sum(sam.lengths)
    
    #select windows
    for bed in bedstream:
        #unload bed coordinate
        chrom, start, end = bed.split('\t')[:3]
        start, end = int(start), int(end)
        #get read count
        c = sam.count(chrom, start, end)
        #check if correct coverage
        cov = c *1.0 / (end-start)
        if cov < covfr * meancov:
            continue
        #report bed
        sys.stdout.write("%s\t%s\t%s\t%.2f\n" % (chrom, start, end, cov))

def main():

    usage  = "%(prog)s [options]"
    parser = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)

    parser.add_argument("-v", "--verbose", default=False, action="store_true")
    parser.add_argument("--version", action="version", version='%(prog)s 0.1')
    parser.add_argument("-a", "--bam", required=True, #default=sys.stdin, type=file,
                        help="input bam")
    parser.add_argument("-b", "--bed", default=sys.stdin, type=file,
                        help="input bed                     [stdin]")
    parser.add_argument("-c", "--covfr", type=float, default=0.75, 
                        help="min fraction of mean coverage [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n" % str(o))
    
    filter_by_coverage(o.bam, o.bed, o.covfr, o.verbose)

if __name__=='__main__':
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!                \n")
    dt = datetime.now()-t0
    sys.stderr.write("## Time elapsed: %s\n" % dt)
