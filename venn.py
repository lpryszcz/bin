#!/usr/bin/env python
desc="""Take list of genes, find overlaps and produce nice Venn.
Dependencies
numpy
pyplot
https://github.com/konstantint/matplotlib-venn
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 7/05/2013
"""

import argparse, os, sys
from datetime import datetime
from Bio      import SeqIO
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles, venn2

def load_columns(files, output, columns, splitname, verbose):
    """Return samples and variants from input files
    Skip #-starting lines
    """
    samples  = []
    variants = []
    for f in files:
        #store name
        name = f.name.split(splitname)[0]
        samples.append(name)
        variants.append(set())
        for l in f:
            l = l[:-1]
            if not l or l.startswith("#"):
                continue
            #unpack columns
            lData = l.split('\t')
            #skip to short columns
            if len(lData)<max(columns):
                sys.stderr.write("Warning: Not enough columns in line: %s\n" % "-".join(lData))
                continue
            #get & store key
            key = "\t".join(lData[c] for c in columns)
            variants[-1].add(key)
    return samples,variants

def variants2venn(samples, variants, output):
    """3-way
    commons need to be: 1, 2, 1-2, 3, 1-3, 2-3, 1-2-3
    or dict
    #for 4-way Venn, need to get another round, but how?
    """
    k = len(variants)
    #get variants stats
    output.write("#Samples\tvariants\n")
    for i,s in enumerate(samples):
        output.write("%s\t%s\n" % (s, len(variants[i])))
    #combine all variants 
    unionset = variants[0]
    for i in range(1,k):
        unionset = unionset.union(variants[i])
    output.write("uniq\t%s\n" % len(unionset))
    venn = {}    
    ##get common pairs
    output.write("#Sample(s)\tcommon\n")
    #get common to all samples
    allcommon = variants[0]
    for i in range(1,k):
        allcommon = allcommon.intersection(variants[i])
    venn['111'] = len(allcommon)
    cname  = " ".join(samples)
    output.write("%s\t%s\n" % (cname, len(allcommon)))
    cnames  = []
    key = ["0","0","0"]
    for i in range(k):
        #get uniq for this sample
        ikey    = key[:]
        ikey[i] = "1"
        sikey   = "".join(ikey)
        #combine all variants
        unionseti = set()
        iis = range(k)
        iis.remove(i)
        for ii in iis:
            unionseti = unionseti.union(variants[ii])
        #get uniq for i
        uniq    = variants[i].difference(unionseti)
        venn[sikey] = len(uniq)
        output.write("%s\t%s\n" % (samples[i], len(uniq)))
        for j in range(i+1,k):
           common  = variants[i].intersection(variants[j]).difference(allcommon)
           #ie 110 or 101
           jkey    = ikey[:] 
           jkey[j] = "1"
           sjkey   = "".join(jkey)
           venn[sjkey] = len(common)
           #get names
           cname = " ".join((samples[i], samples[j]))
           output.write("%s\t%s\n" % (cname, len(common)))                
    return venn

def commons2plot(samples, csizes, output, title, dpi, format, verbose):
    """
    #https://github.com/konstantint/matplotlib-venn
    """
    #get sizes
    plt.figure(figsize=(6,6))
    v = venn3(subsets=csizes, set_labels=samples)
    c = venn3_circles(subsets=csizes, linewidth=0.1) # ,linestyle='dashed'
    plt.title(title)
    #show Venn if not outfile provided
    if output.name=="<stdout>":
        plt.show()
    else:
        fpath = "%s.%s" % (output.name, format)
        plt.savefig(fpath, dpi=dpi, facecolor='w', edgecolor='w',
          orientation='landscape', format=format, transparent=False)
    
def venn(files, output, columns, splitname, title, dpi, format, verbose):
    """
    """
    #load data
    samples,variants = load_columns(files, output, columns, splitname, verbose)

    #get overlaps
    venn = variants2venn(samples, variants, output)
    
    #save venn
    #csizes = (len(c) for c in commons)
    commons2plot(samples, venn, output, title, dpi, format, verbose)

def main():

    usage  = "usage: %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="input",     nargs=3, type=file, 
                        help="input files         [%(default)s]") 
    parser.add_argument("-o", dest="output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream       [stdout]")
    parser.add_argument("-c", dest="columns",   default=[0], type=int, nargs="+",
                        help="data columns        [%(default)s]")
    parser.add_argument("-s", dest="splitname", default=".", 
                        help="split file name     [%(default)s]")
    parser.add_argument("-t", dest="title",     default="Sample Venn", 
                        help="figure title        [%(default)s]")
    parser.add_argument("--dpi", dest="dpi",      default=300, type=int, 
                        help="figure dpi        [%(default)s]" )
    parser.add_argument("--format", dest="format",default="png",
                        choices=['emf', 'eps', 'jpeg', 'jpg', 'pdf', 'png', \
                        'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff'], 
                        help="figure format     [%(default)s]" )
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Args: %s\n" % str(o) )

    venn(o.input, o.output, o.columns, o.splitname, o.title, o.dpi, o.format, o.verbose)

if __name__=='__main__': 
    t0=datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt=datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
