#!/usr/bin/env python3
desc="""Plot scatter plots from many files. First column has to be name of the feature.

TBD
- spaerman
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
from scipy.stats import spearmanr, pearsonr
 
def plot_scatter(fnames, out, col, names, title, log, vmax, vmin, norm, verbose=1):
    """   """
    if not names:
        names = fnames
    if verbose:
        sys.stderr.write( "Parsing data...\n" )
    # get feature positions
    f2i = {l[:-1].split('\t')[0]: i for i, l in enumerate(open(fnames[0], 'r'))}; print(len(f2i))
    data = np.zeros((len(fnames), len(f2i)))
    for i, fn in enumerate(fnames):
        for l in open(fn, 'r'):
            if l.startswith('#'): continue
            ldata = l[:-1].split('\t')
            f = ldata[0]
            if f not in f2i: continue
            try:
                if col>=len(ldata) or not ldata[col]:
                    continue
                v = float(ldata[col])
                if vmin<v<vmax:
                    data[i, f2i[f]] = v
            except:
                sys.stderr.write("[Error] Cannot parse line: %s\n" % ",".join(l[:-1].split('\t')))
    print(data.sum(axis=1)); print(data.shape)
    if norm:
        #data = data* 1e6 / data.sum(axis=1)[:, np.newaxis]
        for i in range(data.shape[0]):
            data[i] *= 1e6 / data.sum(axis=1)[i]
        print(data.sum(axis=1))
    #define number of rows and columns
    ncol = len(fnames)-1
    nrow = len(fnames)-1
    #nrow, ncol = 3, 6
    if verbose:
        sys.stderr.write(" %s columns x %s rows\n"%(ncol, nrow))
    #start figure
    fig, axes = plt.subplots(figsize=(5*ncol, 5*nrow), nrows=nrow, ncols=ncol, sharex=True, sharey=True)
    if title: fig.suptitle(title)
    for i in range(0, len(fnames)-1):
        for j in range(i+1, len(fnames)):
            ax = axes[j-1][i]
            #plot
            x, y = data[i], data[j]; print(np.sum(x>1), np.sum(y>1))
            rho, p = spearmanr(x, y)
            #Prho, Pp = pearsonr(x, y) # Pearson rho=%.4f P=%s
            ax.scatter(x, y, marker='.')
            ax.set_title("%s vs %s\nSpearman rho=%.4f P=%s"%(fnames[i].split(".")[0], fnames[j].split(".")[0], rho, p))#, Prho, Pp))
            #add subplots labels
            if j==nrow:
                ax.set_xlabel(names[i])#, fontsize=30)
            if not i:
                ax.set_ylabel(names[j])#, fontsize=30)
            
    if log:
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlim((1, data.max()))
        ax.set_ylim((1, data.max()))
    #save or show
    if not out:
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
    parser.add_argument("-i", "--input",   nargs='+', help="input files")
    parser.add_argument("-o", "--output",  default=None,  help="output [stdout]")
    parser.add_argument("-c", "--col",     default=0, type=int,
                        help="column to use  [%(default)s]")
    parser.add_argument("-n", "--names", nargs="+", help="column names [%(default)s]")
    parser.add_argument("--log", action="store_true", help="log scale")
    parser.add_argument("-t", "--title", help="figure title")
    parser.add_argument("--max", default=float('inf'), type=float, help="max value [%(default)s]")
    parser.add_argument("--min", default=float('-inf'), type=float, help="min value [%(default)s]")
    parser.add_argument("--norm", action="store_true", help="normalise values [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    plot_scatter(o.input, o.output, o.col, o.names, o.title, o.log, o.max, o.min, o.norm, o.verbose)
              
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
