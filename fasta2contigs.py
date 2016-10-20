#!/usr/bin/env python
desc="""Split FASTA on gaps
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 14/05/2015
"""

import os, re, sys
from datetime import datetime
from optparse import OptionParser
from Bio      import SeqIO

def fasta2contigs(out, handle, lenth, fnformat, verbose):
    """Print BED for every gap above given lenght."""
    #define gap pattern
    gap = re.compile("[Nn]{%s,}" % lenth)
    #store gaps length
    gaps = []
    contigs = set()
    for r in SeqIO.parse(handle, fnformat):
        #iterate gaps
        pe = 0
        for i, m in enumerate(gap.finditer(str(r.seq)), 1):
            s, e = m.span()
            #bed = "%s\t%s\t%s\n" % (r.id, s, e)
            fasta = ">%s.%s\n%s\n"%(r.id, i, str(r.seq)[pe:s])
            out.write(fasta)
            #store info
            pe = e
            gaps.append(e-s)
            contigs.add(r.id)
        # store last bit of given contig
        if pe==0:
            out.write(">%s\n%s\n"%(r.id, str(r.seq)))
        else:
            out.write(">%s.%s\n%s\n"%(r.id, i+1, str(r.seq)[pe:]))
    #print stats
    sys.stderr.write("%s bp in %s gaps in %s contigs.\n" % (sum(gaps), len(gaps), len(contigs)))

def main():
    import argparse
    usage   = "usage: %prog [options]"
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog) 

    parser.add_argument("-v", "--verbose",  default=True, action="store_false")
    parser.add_argument("--version", action='version', version="%(prog)s v1.0")
    parser.add_argument("-i", "--input",   default=sys.stdin, type=file, 
                        help="genome file         [stdin]" )
    parser.add_argument("-o", "--output",   default=sys.stdout, type=argparse.FileType("w"), 
                        help="genome file         [stdout]" )
    parser.add_argument("-f", "--format",   default="fasta", 
                        help="genome file format  [%(default)s]" )
    parser.add_argument("-l", "--gaplength",    default=100, type=int,
                        help="only gaps above     [%(default)s bp]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    #plot
    fasta2contigs(o.output, o.input, o.gaplength, o.format, o.verbose)
  
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)