#!/usr/bin/env python
desc="""2D plot"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 10/05/2013
"""


import argparse, math, os, sys
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

def plot_2d(inputs, output, title, xlab, ylab, xmax, xmin, log, ndivide=20):
    """
    """
    #define number of subplots
    bi = math.sqrt(len(inputs))
    bi = round(bi)
    bj = len(inputs)/bi
    if bj < round(bj):
        bj = round(bj)+1
    else:
        bj = round(bj)+1
    print len(inputs),bi,bj
    #get figure
    plt.figure(figsize=(bj*4, bi*4))
    plt.subplots_adjust(hspace = .3, wspace = .3)    
    #process inputs    
    sys.stderr.write("Loading data...\n")
    for ii,input in enumerate(inputs, 1):
        #load data
        x, y = [], []
        for l in input:
            l = l[:-1]
            if not l or l.startswith('#'):
                continue
            i, c = l.split()[:2]
            i, c = int(i), int(c)
            if xmin <= i <= xmax:
                x.append(i)
                y.append(c/10**3)

        maxy10 = max(y[xmax/ndivide:])
        xi = y.index(maxy10)
        freqk = x[xi]
        gsize = maxy10*freqk/10.0**3
        sys.stderr.write("[%s] Max freq: %s @ k-mer freq: %s\nEstimated genome size: %s Mb\n" %(input.name, maxy10, freqk, gsize))
        plt.subplot(bi,bj,ii)
        plt.plot(x, y, linewidth=2.0)
        #add title and axis labels
        if input.name!="<stdin>":
            plt.title(input.name.split('.')[0])
        elif title:
            plt.title(input.name)
        #plot x-axis label only on bottom plots
        if ii/bi>bj-1:
            plt.xlabel(xlab)
        #plot y-axis label only on left-most plots
        if ii%bj==1:
            plt.ylabel(ylab)
        plt.ylim(0,1.5*maxy10)
        #plt.grid(True)
        #add local max
        plt.annotate("~%.2f Mb\n(%s, %sK)" % (gsize, freqk, maxy10), xy=(freqk*1.01, maxy10*1.01), xytext=(freqk*1.2, maxy10*1.2),arrowprops=dict(facecolor='black', shrink=0.05))
        #plt.text(freqk, maxy10*1.1, 'Genome size: ~%.2f Mb\n(%s, %s)' % (gsize, freqk, maxy10))
    #show plot if not outfile provided    
    if output.name=="<stdout>":
        plt.show()
    else:
        fpath  = output.name #"%s.%s" % (output.name, format)
        format = fpath.split('.')[-1]
        plt.savefig(fpath, dpi=200, facecolor='w', edgecolor='w',\
          orientation='landscape', format=format, transparent=False)
        
def main():
    
    usage   = "%(prog)s [options] -v" 
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
  
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="input",   default=[sys.stdin,], type=file, nargs="+",
                        help="input stream    [stdin]")
    parser.add_argument("-o", dest="output",  default=sys.stdout, type=argparse.FileType("w"), 
                        help="output stream   [stdout]")
    parser.add_argument("-c", dest="col",     default=0, type=int,
                        help="column to use   [%(default)s]")
    parser.add_argument("-t", dest="title",   default="", 
                        help="histogram title [%(default)s]")
    parser.add_argument("-x", dest="xlab",    default="k-mer frequency", 
                        help="x-axis label    [%(default)s]")
    parser.add_argument("-y", dest="ylab",    default="k-mers with this frequency [10e3]", 
                        help="y-axis label    [%(default)s]")
    parser.add_argument("-n", dest="ndivide", default=20, type=int,
                        help="discard 1/n first     [%(default)s]")
    parser.add_argument("--log", dest="log",  default=False, action="store_true",
                        help="log scale       [%(default)s]")
    parser.add_argument("--xmax", dest="xmax", default=100, type=int,
                        help="max x value     [%(default)s]")
    parser.add_argument("--xmin", dest="xmin", default=0, type=int,
                        help="min x value     [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    plot_2d(o.input, o.output, o.title, o.xlab, o.ylab, o.xmax, o.xmin, o.log, o.ndivide)
              
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
