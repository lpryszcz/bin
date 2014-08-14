#!/usr/bin/env python
"""
Parse BED and print stats.

Author:
l.p.pryszcz@gmail.com

Dublin, 4/09/2012
"""

import os, sys
import numpy as np
from optparse import OptionParser,OptionGroup
from datetime import datetime
from genome_annotation import get_contig2coverage

import numpy as np
import matplotlib.pyplot as plt

def print_stats( lengths,addon="" ):
    """
    """
    info = "%.2f kb in %s fragments%s (median: %s bp; mean: %.3f kb +-%.3f )"
    print info%(sum(lengths)/10.0**3, len(lengths), addon, \
                int(np.median(lengths)), np.mean(lengths)/10.0**3,
                np.std(lengths)/10.0**3)

def bed2stats( handle,simple,verbose ):
    """Parse BED and print stats."""
    r2l = {}
    for l in handle:
        ref,start,end = l.split()[:3]
        start,end = int(start),int(end)
        r2l[ "%s:%s-%s" % (ref,start,end) ] = end-start
    #print summary
    lengths = [ x for x in r2l.itervalues() ]
    if not simple:    
        print_stats( lengths )
    else:
        print "%s\t%s\t%s" % ( handle.name,sum(lengths),len(r2l) )
        return
    
    longest = ""
    r2l_sorted = sorted( r2l.iteritems(),key=lambda x: x[1], reverse=True )
    for r,l in r2l_sorted[:20]:
        longest += "\t\t%s\t%s\n" % ( l,r )
    longest += "\t\t...\n\t\t%s\t%s" % ( r2l_sorted[-1][1],r2l_sorted[-1][0] )
    if not simple:
        print longest

    lengths1kb = []
    for x in r2l.itervalues():
        if x>=1000:
            lengths1kb.append( x )
    if not simple:
        print_stats( lengths1kb," >=1kb" )    

def main():
    usage  = "usage: bedtools merge -d 1 -i bed | %prog [options]"
    desc   = "Parse BED and print stats."
    epilog = ""
    parser = OptionParser( usage=usage,version="%prog 1.0",description=desc,epilog=epilog ) 

    parser.add_option("-v", dest="verbose", default=False, action="store_true" )
    #parser.add_option("-a", dest="bam",  
    #                  help="BAM file               [mandatory]")
    parser.add_option("-b", dest="bed", default="",
                      help="BED file               [stdin]")   
    #parser.add_option("-c", dest="cov_fract", default=0.75, type=float,
    #                  help="frac of mean coverage  [%default]")
    parser.add_option("-s", dest="simple", default=False, action="store_true",
                      help="simple output           [%default]")       
    
    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nArgs: %s\n" % ( o,args ) )

    fnames = args
    if o.bed:
        fnames = [ o.bed, ]
    if o.simple:
        print "#sample\tsum\toccurencies"
    for fn in fnames:
        if not fn:
            handle = sys.stdin
        elif not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )
        else:
            handle = open( fn )
    
        bed2stats( handle,o.simple,o.verbose )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  #sys.stderr.write( "#Time elapsed: %s\n" % dt )
