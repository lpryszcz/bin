#!/usr/bin/env python
desc="""Plot histogram
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 31/07/2013
"""

import argparse, os, sys
from datetime import datetime

def load_idxstats(input):
    """Return chr2len and chr2reads."""
    chr2len   = {}
    chr2reads = {}
    for l in input:
        chrom, size, reads, unaligned = l.split('\t')
        if chrom=="*":
            continue
        chr2len[chrom]   = int(size)
        chr2reads[chrom] = int(reads)
    return chr2len, chr2reads

def chr2copynumber(chrom, chr2len, chr2reads, ploidy):
    """Return copy number of give chromosome."""
    #genome size and total reads
    gsize  = sum(chr2len.itervalues())
    treads = sum(chr2reads.itervalues())
    #expected read count
    ereads = chr2len[chrom] * treads * 1.0 / gsize
    #print ereads, chr2reads[chrom]
    #copy number
    return ploidy * chr2reads[chrom] / ereads
    
def idxstats2copynumber(input, ploidy, chromosomes, verbose):
    """Report chromosome copy number."""
    chr2len, chr2reads = load_idxstats(input)

    #select all chromosomes if none specified
    if not chromosomes:
        chromosomes = sorted(chr2len.keys(), key=lambda x: chr2len[x], reverse=1)
    #process all chromosomes
    for chrom in chromosomes:
        if chrom not in chr2len:
            sys.stderr.write("[Error] %s not in chromosomes: %s\n" % (chrom, ",".join(chr2len.keys()[:10])))
        chromcn = chr2copynumber(chrom, chr2len, chr2reads, ploidy)
        sys.stdout.write("%s\t%5.2f\n" % (chrom, chromcn))

def main():
    
    usage   = "%(prog)s [options] -v" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v","--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="input",   default=sys.stdin, type=file,
                        help="input           [stdin]")
    parser.add_argument("--ploidy",           default=2, type=int,
                        help="ploidy    [%(default)s]")
    parser.add_argument("-c", "--chroms", nargs="*", 
                        help="chromosome(s) to analyse [all]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    idxstats2copynumber(o.input, o.ploidy, o.chroms, o.verbose)
              
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    #sys.stderr.write( "#Time elapsed: %s\n" % dt )
