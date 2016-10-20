#!/usr/bin/env python
desc="""Parse BED and print stats.
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Dublin, 4/09/2012
"""

import os, sys
import numpy as np
#from optparse import OptionParser,OptionGroup
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
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-b", "--bed", nargs="+", default=[sys.stdin], type=file, 
                        help="BED file        [stdin]")   
    #parser.add_option("-c", dest="cov_fract", default=0.75, type=float,
    #                  help="frac of mean coverage  [%default]")
    parser.add_argument("-s", dest="simple", default=False, action="store_true",
                        help="simple output   [%(default)s]")       
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    if o.simple:
        print "#sample\tsum\toccurencies"
    for handle in o.bed:
        bed2stats(handle, o.simple, o.verbose)

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
