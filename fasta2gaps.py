#!/usr/bin/env python2
desc="""Report BED of gaps for given genome
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 3/10/2012
"""

import os, re, sys
from datetime import datetime
from optparse import OptionParser
from Bio      import SeqIO

def genome2gaps(fn, lenth, fnformat, verbose):
    """Print BED for every gap above given lenght."""
    #define gap pattern
    gap = re.compile("[Nn]{%s,}" % lenth)
    #store gaps length
    gaps = []
    contigs = set()
    for r in SeqIO.parse(open(fn), fnformat):
        #iterate gaps
        for m in gap.finditer(str(r.seq)):
            s,e = m.span()
            bed = "%s\t%s\t%s\n" % (r.id, s, e)
            sys.stdout.write(bed)
            #store info
            gaps.append(e-s)
            contigs.add(r.id)
    #print stats
    sys.stderr.write("%s bp in %s gaps in %s contigs.\n" % (sum(gaps), len(gaps), len(contigs)))

def main():
    usage   = "usage: %prog [options]"
    version = "%prog 1.0"
    parser  = OptionParser(usage=usage, version=version, description=desc, epilog=epilog) #allow_interspersed_args=True

    parser.add_option("-v", dest="verbose",  default=True, action="store_false")
    parser.add_option("-i", dest="genome",   default="", 
                      help="genome file         [%default]" )
    parser.add_option("-f", dest="format",   default="fasta", 
                      help="genome file format  [%default]" )
    parser.add_option("-l", dest="length",    default=100, type=int,
                      help="only gaps above     [%default bp]")
    
    o,args = parser.parse_args()
    if o.verbose: 
        sys.stderr.write( "Options: %s\nArgs: %s\n" % (str(o),", ".join(args)))

    if not o.genome:
        parser.error("Specify genome file!")
    if not os.path.isfile( o.genome ):
        parser.error("No such file: %s!" % o.genome)

    if o.length<1:
        parser.error("Gap length has to be an positive integer!")
      
    #plot
    genome2gaps(o.genome, o.length, o.format, o.verbose)
  
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)