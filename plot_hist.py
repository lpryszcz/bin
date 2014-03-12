#!/usr/bin/env python
desc="""Plot histogram
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 29/11/2012
"""

import argparse, os, sys
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

def histplot(out, x, handle, categories, title, xlab, ylab, log):
    """
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if log:
    	ax.set_yscale('log')
    # the histogram of the data
    n, bins, patches = ax.hist(x, categories, normed=0, facecolor='green', alpha=0.75)

    # hist uses np.histogram under the hood to create 'n' and 'bins'.
    # np.histogram returns the bin edges, so there will be 50 probability
    # density values in n, 51 bin edges in bins and 50 patches.  To get
    # everything lined up, we'll compute the bin centers
    #bincenters = 0.5*(bins[1:]+bins[:-1])
    # add a 'best fit' line for the normal PDF
    #y = mlab.normpdf( bincenters, mu, sigma)
    #l = ax.plot(bincenters, y, 'r--', linewidth=1)

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title+'\nmean=%.2f; stdev=%.2f' % (np.mean(x),np.std(x)))
    #ax.set_title('mean=%.2f,stdev=%.2f' % (np.mean(x),np.std(x)))
    #ax.set_xlim(40, 160)
    #ax.set_ylim(0, 0.03)
    #plt.text(60, .025, r'mean=%.2f,stdev=%.2f' % (np.mean(x),np.std(x)))
    ax.grid(True)
    if type(out) is file and out.name=='<stdout>':
    	plt.show()
    else:
        fpath = out #handle.name+".png"
        fformat = fpath.split('.')[-1]
        plt.savefig(fpath, orientation='landscape', format=fformat, transparent=False)
        print "Figure written to: %s" % fpath
 

def plot_hist(handle, out, col, bins, title, xlab, ylab, log, vmax, vmin, verbose):
    """
    """
    if verbose:
        sys.stderr.write( "Parsing data...\n" )
    x = []
    for l in handle:
        l = l.strip()
        if not l:
            continue
        try:
            v=float(l.split('\t')[col])
            if vmin<v<vmax:
        	x.append(v)
        except:
            sys.stderr.write("[Error] Cannot parse line: %s\n" % ", ".join(l.split('\t')))
    if verbose:
        sys.stderr.write( " %s values loaded.\n" % len(x) )
    #print( x[:10],min(x),max(x) )

    histplot(out, x, handle, bins, title, xlab, ylab, log)
    
def main():
    
    usage   = "%(prog)s [options] -v" 
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
  
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="input",   default=sys.stdin, type=file,
                        help="input           [stdin]")
    parser.add_argument("-o", dest="output",  default=sys.stdout, 
                        help="input           [stdout]")
    parser.add_argument("-b", dest="bins",    default=100, type=int,
                        help="number of bins  [%(default)s]")
    parser.add_argument("-c", dest="col",     default=0, type=int,
                        help="column to use   [%(default)s]")
    parser.add_argument("-t", dest="title",   default="Histogram", 
                        help="histogram title [%(default)s]")
    parser.add_argument("-x", dest="xlab",    default="", 
                        help="x-axis label    [%(default)s]")
    parser.add_argument("-y", dest="ylab",    default="", 
                        help="y-axis label    [%(default)s]")
    parser.add_argument("--log", dest="log",  default=False, action="store_true",
                        help="log scale       [%(default)s]")
    parser.add_argument("--max", dest="max",  default=float('inf'), type=float,
                        help="max value       [%(default)s]")
    parser.add_argument("--min", dest="min",  default=float('-inf'), type=float,
                        help="min value       [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    plot_hist(o.input, o.output, o.col, o.bins, o.title, o.xlab, o.ylab, o.log, \
              o.max, o.min, o.verbose)
              
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
