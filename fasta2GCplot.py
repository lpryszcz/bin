#!/usr/bin/env python
desc="""Report GC content.
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 12/11/2013
"""

import argparse, gzip, os, sys
from datetime import datetime
from Bio import SeqIO
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import numpy as np

def fasta2gc(fasta, verbose):
    """Return sequences GC and length"""
    handle = fasta
    if fasta.name.endswith('.gz'):
        handle = gzip.open(fasta)
    #id2gc = {}
    gc, lengths = [], []
    if verbose:
        sys.stderr.write("Processing %s ...\n"%fasta.name)
    for r in SeqIO.parse(handle, 'fasta'):
        if verbose:
            sys.stderr.write(" %s [%s bp]: %.2f%s\n" %(r.id, len(r), GC(r.seq), '%'))
        gc.append(GC(r.seq))
        lengths.append(len(r))
    return gc, lengths #id2gc

def GCplot(fastas, ofpath, title, xlab, ylab, log, verbose):
    """ """
    #init plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if log:
        ax.set_yscale('log')    
        
    #get GC
    colors = ('bo','ro','go','mo','yo') #('bo','r>','gv','m^','y<')
    for handle, c in zip(fastas, colors):
        gc, lengths = fasta2gc(handle, verbose)
        ax.plot(gc, lengths, c, alpha=0.75, label=handle.name)
    #plot
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.grid(True)
    ax.legend()
    if not ofpath:
    	plt.show()
    else:
        fformat = ofpath.split('.')[-1]
        plt.savefig(ofpath, orientation='landscape', format=fformat, transparent=False)
        print "Figure written to: %s" % ofpath
    
def main():
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",   default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-f", "--fastas",     nargs="+", type=file, 
                        help="fasta stream(s)")
    parser.add_argument("-o", "--output",     default=0, 
                        help="output   [X11]")
    parser.add_argument("-t", dest="title",   default="Plot", 
                        help="figure title    [%(default)s]")
    parser.add_argument("-x", "--xlab",        default="GC[%]", 
                        help="x-axis label    [%(default)s]")
    parser.add_argument("-y", "--ylab",        default="length", 
                        help="y-axis label    [%(default)s]")
    parser.add_argument("--log",               default=False, action="store_true",
                        help="log scale")
    '''parser.add_argument("--max", dest="max",  default=float('inf'), type=float,
                        help="max value       [%(default)s]")
    parser.add_argument("--min", dest="min",  default=float('-inf'), type=float,
                        help="min value       [%(default)s]")'''
                         
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n" % str(o))

    GCplot(o.fastas, o.output, o.title, o.xlab, o.ylab, o.log, o.verbose)
	
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
