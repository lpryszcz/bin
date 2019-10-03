#!/usr/bin/env python3
desc="""Plot histogram from many files
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 29/11/2012
"""

import argparse, os, sys
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from io import IOBase # type file

def histplot(ax, x, bins, title, xlab, ylab, xlog, ylog, \
             normed=False, cumulative=False, alpha=0.75):
    """
    """
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    if xlog:
    	ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log', nonposy='clip')
    # the histogram of the data
    n, bins, patches = ax.hist(x, bins, density=normed, alpha=alpha,
                               cumulative=cumulative, label=title)

    # hist uses np.histogram under the hood to create 'n' and 'bins'.
    # np.histogram returns the bin edges, so there will be 50 probability
    # density values in n, 51 bin edges in bins and 50 patches.  To get
    # everything lined up, we'll compute the bin centers
    #bincenters = 0.5*(bins[1:]+bins[:-1])
    # add a 'best fit' line for the normal PDF
    #y = mlab.normpdf( bincenters, mu, sigma)
    #l = ax.plot(bincenters, y, 'r--', linewidth=1)
    #title+'\nmean=%.2f; stdev=%.2f' % (np.mean(x),np.std(x)), \
    ax.set_title(title)
    #ax.set_title('mean=%.2f,stdev=%.2f' % (np.mean(x),np.std(x)))
    #ax.set_xlim(40, 160)
    #ax.set_ylim(0, 0.03)
    #plt.text(60, .025, r'mean=%.2f,stdev=%.2f' % (np.mean(x),np.std(x)))
    ax.grid(True)
    return bins
 
def plot_hist(fnames, out, col, names, bins, title, xlab, ylab, xlog, ylog, \
              vmax, vmin, collapse, normed, cumulative, alpha, legendLoc, colors,\
              figsize=(15,15), verbose=1, dlimit=1):
    """
    """
    print(col, fnames)
    if not names:
        names=fnames
    if verbose:
        sys.stderr.write( "Parsing data...\n" )
    x = cols = [[] for fn in fnames]
    for i, fn in enumerate(fnames):
        for l in open(fn, 'r'):
            if l.startswith('#'): continue
            try:
                ldata = l[:-1].split('\t')
                if col>=len(ldata) or not ldata[col]:
                    continue
                v=float(ldata[col])
                if vmin<v<vmax:
                    x[i].append(v)
            except:
                sys.stderr.write("[Error] Cannot parse line: %s\n" % ",".join(l[:-1].split('\t')))
    if verbose:
        sys.stderr.write( " %s values loaded.\n" % len(x) )
    #define number of rows and columns
    ncol = nrow = int(np.sqrt(len(cols)))
    if ncol*nrow < len(cols):
        ncol += 1
    if ncol*nrow < len(cols):
        nrow += 1
    #nrow, ncol = 3, 6
    if verbose:
        sys.stderr.write(" %s columns x %s rows\n"%(ncol, nrow))
    #start figure
    fig = plt.figure(figsize=figsize) # set size
    #http://matplotlib.org/users/customizing.html
    #mpl.rcParams['figure.subplot.wspace'] = 0.3
    mpl.rcParams['figure.subplot.hspace'] = 0.5
    mpl.rcParams['axes.titlesize'] = 8
    mpl.rcParams['axes.labelsize'] = 6
    mpl.rcParams['xtick.labelsize'] = 5
    mpl.rcParams['ytick.labelsize'] = 5
    #add subplots
    #plt.rc('axes', color_cycle=colors) #['c', 'm', 'y', 'k']
    #mpl.rcParams['axes.color_cycle'] = colors
    for i, data in enumerate(x):
        if collapse:
            if i==0:
                mpl.rcParams['axes.titlesize'] = 18
                mpl.rcParams['axes.labelsize'] = 10
                mpl.rcParams['xtick.labelsize'] = 9
                mpl.rcParams['ytick.labelsize'] = 9
                ax = fig.add_subplot(111)
        else:
            ax = fig.add_subplot(nrow, ncol, i+1)
        #get name
        name = ""
        if names:
            name = names[i]
        if len(data)<dlimit:
            sys.stderr.write("[WARNING] Only %s data points for: %s\n"%(len(data), name))
            continue
        #plot
        bins = histplot(ax, data, bins, name, xlab, ylab, xlog, ylog, normed, cumulative, alpha)
        #add title
        if len(cols)==1 and title or collapse:
            ax.set_title(title)
        #add subplots labels
        if i+1>=len(cols)-nrow:
            ax.set_xlabel(xlab)#, fontsize=30)
        if not i%ncol:
            ax.set_ylabel(ylab)#, fontsize=30)
    #plot legend only if collapsed
    if collapse:
        ax.legend(loc=legendLoc)
    #save or show
    if isinstance(out, IOBase) and out.name=='<stdout>':
    	plt.show()
    else:
        fpath = out #handle.name+".png"
        fformat = fpath.split('.')[-1] 
        plt.savefig(fpath, dpi=300, format=fformat, orientation='landscape', transparent=False)
        sys.stderr.write("Figure written to: %s\n" % fpath)
    
    
def main():
    
    usage   = "%(prog)s [options] -v" 
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
  
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", "--input",   nargs='+', #type=argparse.FileType("r"),
                        help="input files")
    parser.add_argument("-o", "--output",  default=sys.stdout, 
                        help="input           [stdout]")
    parser.add_argument("-b", "--bins",    default=100, type=int,
                        help="number of bins  [%(default)s]")
    parser.add_argument("-c", "--col",     default=0, type=int,
                        help="column to use  [%(default)s]")
    parser.add_argument("-n", "--names",   default=None, nargs="+", 
                        help="column names    [%(default)s]")
    parser.add_argument("-t", "--title",   default="Histogram", 
                        help="histogram title [%(default)s]")
    parser.add_argument("-x", "--xlab",    default="", 
                        help="x-axis label    [%(default)s]")
    parser.add_argument("-y", "--ylab",    default="", 
                        help="y-axis label    [%(default)s]")
    parser.add_argument("--ylog", "--log", default=False, action="store_true",
                        help="log Y scale")
    parser.add_argument("--xlog",          default=False, action="store_true",
                        help="log X scale")
    parser.add_argument("--max",           default=float('inf'), type=float,
                        help="max value       [%(default)s]")
    parser.add_argument("--min",           default=float('-inf'), type=float,
                        help="min value       [%(default)s]")
    parser.add_argument("--collapse",      default=False, action="store_true",
                        help="collapse into single subplot")
    parser.add_argument("--cumulative",  default=0, choices=[0, 1, -1], type=int, 
                        help="cumulative histogram [%(default)s]")
    parser.add_argument("--normed",        default=0, choices=(0, 1), type=int,
                        help="normalise values")
    parser.add_argument("--alpha",         default=0.75, type=float,
                        help="plot alpha      [%(default)s]")
    parser.add_argument("--figsize", default=(15, 15), nargs=2, type=int, 
                        help="figure dimensions [%(default)s]")
    parser.add_argument("--colors",        nargs="+", default=['b', 'grey', 'r', 'y', 'g'],
                        help="plot alpha      [%(default)s]")    
    parser.add_argument("--legendLoc",     default=1, choices=(1, 2, 3, 4), type=int,
                        help="legend location (1=top right, 2=top left, 3=bottom left, 4=bottom right")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    plot_hist(o.input, o.output, o.col, o.names, o.bins, o.title, o.xlab, o.ylab, \
              o.xlog, o.ylog, o.max, o.min, o.collapse, o.normed, o.cumulative, o.alpha, o.legendLoc, \
              o.colors, o.figsize, o.verbose)
              
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
