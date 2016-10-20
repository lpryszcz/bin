#!/usr/bin/env python
desc="""Plot histogram
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

import plotly.plotly as py
import plotly.graph_objs as go

data = [
    go.Heatmap(
        z=[[1, 20, 30, 50, 1], [20, 1, 60, 80, 30], [30, 60, 1, -10, 20]],
        x=['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'],
        y=['Morning', 'Afternoon', 'Evening']
    )
]
plot_url = py.plot(data, filename='labelled-heatmap')

def histplot(ax, x, handle, bins, title, xlab, ylab, xlog, ylog, \
             normed=False, alpha=0.75):
    """
    """
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    if xlog:
    	ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log', nonposy='clip')
    # the histogram of the data
    n, bins, patches = ax.hist(x, bins, normed=normed, alpha=alpha, label=title)

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
 
def plot_hist(handle, out, cols, names, bins, title, xlab, ylab, xlog, ylog, \
              vmax, vmin, collapse, normed, alpha, legendLoc, colors,\
              verbose, dlimit=1):
    """
    """
    if verbose:
        sys.stderr.write( "Parsing data...\n" )
    x = [[] for i in range(len(cols))]
    for l in handle:
        try:
            ldata = l[:-1].split('\t')
            for i, col in enumerate(cols):
                if col>=len(ldata) or not ldata[col]:
                    continue
                v=float(ldata[col])
                if vmin<v<vmax:
                    x[i].append(v)
        except:
            sys.stderr.write("[Error] Cannot parse line: %s\n" % ",".join(l.split('\t')))
    if verbose:
        sys.stderr.write( " %s values loaded.\n" % len(x) )
    #define number of rows and columns
    ncol = nrow = int(np.sqrt(len(cols)))
    if ncol*nrow < len(cols):
        ncol += 1
    if ncol*nrow < len(cols):
        nrow += 1
    nrow, ncol = 3, 6
    if verbose:
        sys.stderr.write(" %s columns x %s rows\n"%(ncol, nrow))
    #start figure
    fig = plt.figure()
    #http://matplotlib.org/users/customizing.html
    #mpl.rcParams['figure.subplot.wspace'] = 0.3
    mpl.rcParams['figure.subplot.hspace'] = 0.5
    mpl.rcParams['axes.titlesize'] = 8
    mpl.rcParams['axes.labelsize'] = 6
    mpl.rcParams['xtick.labelsize'] = 5
    mpl.rcParams['ytick.labelsize'] = 5
    #add subplots
    plt.rc('axes', color_cycle=colors) #['c', 'm', 'y', 'k']
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
        bins = histplot(ax, data, handle, bins, name, xlab, ylab, xlog, ylog, normed)
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
    if type(out) is file and out.name=='<stdout>':
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
    parser.add_argument("-i", "--input",   default=sys.stdin, type=file,
                        help="input           [stdin]")
    parser.add_argument("-o", "--output",  default=sys.stdout, 
                        help="input           [stdout]")
    parser.add_argument("-b", "--bins",    default=100, type=int,
                        help="number of bins  [%(default)s]")
    parser.add_argument("-c", "--col",     default=[0], nargs="+", type=int,
                        help="columns to use  [%(default)s]")
    parser.add_argument("-n", "--names",   default="", nargs="+", 
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
    parser.add_argument("--normed",        default=0, choices=(0, 1), type=int,
                        help="normalise values")
    parser.add_argument("--alpha",         default=0.75, type=float,
                        help="plot alpha      [%(default)s]")
    parser.add_argument("--colors",        nargs="+", default=['b', 'grey', 'r', 'y', 'g'],
                        help="plot alpha      [%(default)s]")    
    parser.add_argument("--legendLoc",     default=1, choices=(1, 2, 3, 4), type=int,
                        help="legend location (1=top right, 2=top left, 3=bottom left, 4=bottom right")
    

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    plot_hist(o.input, o.output, o.col, o.names, o.bins, o.title, o.xlab, o.ylab, \
              o.xlog, o.ylog, o.max, o.min, o.collapse, o.normed, o.alpha, o.legendLoc, \
              o.colors, o.verbose)
              
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
