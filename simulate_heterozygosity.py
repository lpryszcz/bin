#!/usr/bin/env python
desc="""Simulate R heterozygous genome(s) and report their statistics. 
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow, 5/03/2015
"""

import argparse, os, random, sys
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

def get_random_sample(population, divergence):
    """Return random genome positions."""
    # get sample size
    sampleSize = int(round(0.045*population))
    # return random sample
    return random.sample(xrange(population), sampleSize)

def get_distances(sample, pp=None):
    """Return distances"""
    sample.sort()
    d = []
    for p in sample:
        if pp:
            d.append(p-pp)
        pp = p
    return d

def get_loh(distances, minl=100):
    """Return LOH for distances"""
    loh = [0]
    for d in distances:
        if d>=minl:
            loh[-1]+=d
        else:
            if loh[-1]:
                loh.append(0)
    return loh
    
def get_stats(gsize, divergence, repeats, verbose):
    """Report statistics of random heterozygous genome(s)"""
    xlab = "LOH size"
    ylab = "Frequency"
    info = '%s - %s bp\n>=100bp:\n  %s bp\n  %s LOHs\n  mean %.1f\n  median: %i'
    # define number of rows and columns
    ncol = nrow = int(np.sqrt(repeats))
    if ncol*nrow < repeats:
        ncol += 1
    if ncol*nrow < repeats:
        nrow += 1
    #start figure
    fig = plt.figure()    
    for i in range(repeats):
        # get random genome positions
        sample = get_random_sample(gsize, divergence)
        # call LOH
        d = get_distances(sample)
        d100 = get_loh(d, 100) #filter(lambda x: x>=100, d)
        # plot hist
        ax = fig.add_subplot(nrow, ncol, i+1)
        ax.hist(d, 50, normed=1)
        # add info
        plt.text(50, .01, info%(min(d), max(d), sum(d100), len(d100),
                                  np.mean(d100), np.median(d100)))
        print i+1, info%(min(d), max(d), sum(d100), len(d100),
                         np.mean(d100), np.median(d100))
        #add subplots labels
        if i+1>=repeats-nrow:
            ax.set_xlabel(xlab)
        if not i%ncol:
            ax.set_ylabel(ylab)
        #set the same limits
        ax.set_xlim(0, 150)
        ax.set_ylim(0, 0.05)
            
    plt.show()
    
def main():

    usage  = "%(prog)s [options]"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-g", "--gsize", type=int, default=13393615, 
                        help="genome size in bp [%(default)s]")
    parser.add_argument("-d", "--divergence", type=float, default=0.045,
                        help="divergence between haplotypes [%(default)s]")
    parser.add_argument("-r", "--repeats", type=int, default=1,
                        help="repeat r times [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    get_stats(o.gsize, o.divergence, o.repeats, o.verbose)
        
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
    