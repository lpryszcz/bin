#!/usr/bin/env python
"""
Parse BED and print only passing count criteria.

Author:
l.p.pryszcz@gmail.com

Dublin, 4/09/2012
"""

import os, sys
import numpy as np
from optparse import OptionParser,OptionGroup
from datetime import datetime
from genome_annotation import get_contig2coverage

def bed2count_filter( handle,minC,maxC,minL,verbose ):
    """Parse BED and print stats."""
    i = skipped = 0
    for l in handle:
        i += 1
        c = 0
        ref,start,end = l.split()[:3]
        if len(l.split())>3:
            c = float(l.split()[3])
        start,end = int(start),int(end)
        if minC <= c <= maxC and end-start>=minL:
            sys.stdout.write( l )
        else:
            skipped += 1
    if verbose:
        sys.stderr.write( "%s entries of which %s skipped [%.2f%s].\n" % (i,skipped,skipped*100.0/i,'%') )
        
def main():
    usage  = "usage: bedtools intersect -c -a bed1 -b bed2 | %prog [options]"
    desc   = "Parse BED and print entries passing count criteria."
    epilog = ""
    parser = OptionParser( usage=usage,version="%prog 1.0",description=desc,epilog=epilog ) 

    #parser.add_option("-a", dest="bam",  
    #                  help="BAM file               [mandatory]")
    parser.add_option("-b", dest="bed",
                      help="BED file           [stdin]")   
    parser.add_option("-c", dest="minC", default='-inf', type=float,
                      help="min count allowed  [%default]")
    parser.add_option("-d", dest="maxC", default='+inf', type=float,
                      help="max count allowed  [%default]")
    parser.add_option("-l", dest="minL", default=0, type=int,
                      help="min entry length allowed  [%default]")    
    parser.add_option("-v", dest="verbose", default=False, action="store_true" )
    
    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nArgs: %s\n" % ( o,args ) )

    for fn in [ o.bed, ]:
        if not fn:
            handle = sys.stdin
        elif not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )
        else:
            handle = open( o.bed )
    
    bed2count_filter( handle,o.minC,o.maxC,o.minL,o.verbose )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  #sys.stderr.write( "#Time elapsed: %s\n" % dt )
