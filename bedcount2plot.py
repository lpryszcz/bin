#!/usr/bin/env python2
desc="""Generate plot of raw counts from BED files. 

First generate cov.bed files is:
for f in *.bam; do if [ ! -s $f.gff.bed ]; then echo `date` $f; bedtools coverage -abam $f -b yeast_chromosomes.gff > $f.gff.bed; fi; done

Note: outname has to finish with one of the following: png, pdf, ps, eps and svg
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 14/03/2013
"""

import argparse, os, sys
from datetime import datetime
#from math     import log
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr

colors = ('b','r','g','y')

def load_counts(handle, chr2counts, countCol):
    """ """
    #populate empty list to each chr
    for chrname in chr2counts:
        chr2counts[chrname].append([])
    #populate counts
    for line in handle:
        if not line or line.startswith('#'):
            continue
        lData = line.split('\t')
        #skip no genes
        if lData[2]!="gene":
            continue
        #check if count present
        if len(lData)<countCol:
            sys.exit("ERROR: No count column (%s) in: %s" % (lData,handle.name))
        chrname = lData[0]
        #add chr info
        if chrname not in chr2counts:
            chr2counts[chrname]=[[]]
        count = int(lData[countCol])
        #add count
        chr2counts[chrname][-1].append(count)
    return chr2counts

def normalise_sample(chr2counts, ratios):
    """Normalise by total number of reads"""
    #normalise xbed counts
    for chrname in chr2counts:
        for i,ratio in enumerate(ratios):
            #print chrname,i,ratio
            #print chr2counts[chrname][i][:10]
            chr2counts[chrname][i] = [x*ratio for x in chr2counts[chrname][i]]
            #print chr2counts[chrname][i][:10]
    return chr2counts
    
def bedcount2plot(xbeds, ybeds, outname, logscale, normalise, countCol, dpi, verbose, marker="."):
    """
    , marker is one pixel
    """
    #load bed files
    if verbose:
        sys.stderr.write("Loading gene counts from BEDs...\n")
    x2counts, y2counts = {}, {}
    for xbed, ybed in zip(xbeds, ybeds):
        #get gene counts
        x2counts = load_counts(xbed, x2counts, countCol)
        y2counts = load_counts(ybed, y2counts, countCol)
        
    #normalise y by xbed total read count
    if verbose:
        sys.stderr.write("Normalising...\n")
    if normalise:
        #ncount = sum( sum(xlist[0]) for chrname,xlist in x2counts.iteritems() )
        #get samples ratios
        xratios = []
        yratios = [] 
        #normalise all counts
        for i in range(len(zip(xbeds, ybeds))):
            #append x sample ratio
            xscount = sum(sum(xlist[i]) for chrname,xlist in x2counts.iteritems())
            #xratios.append( 1.0*ncount/xscount )
            #append y sample ratio
            yscount = sum(sum(ylist[i]) for chrname,ylist in y2counts.iteritems())
            yratios.append( 1.0*xscount/yscount )
            print xscount,yscount
        #normalise counts
        #x2counts = normalise_sample(x2counts, xratios)
        y2counts = normalise_sample(y2counts, yratios)
        #print y2counts
        print yratios
    #plot
    if verbose:
        sys.stderr.write("Plotting...\n")
    i = 0
    fig = plt.figure()
    X,Y = [],[]      
    fig.subplots_adjust(hspace=.3)
    chrnames = ("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI") #sorted(x2counts.keys()):
    for chrname in chrnames:
        sys.stderr.write("  %s      \r" % chrname)
        i += 1
        #new subplot
        plt.subplot(4,4,i)
        title = chrname+" r: "
        for c,x,y in zip(colors,x2counts[chrname],y2counts[chrname]):
            X += x
            Y += y
            #print chrname, len(x),len(y)
            # red dashes, blue squares and green triangles
            plt.plot(x, y, c+marker, markersize=0.10)
            #get correlation
            pearR,pval = pearsonr(x,y)
            title += "%.5f "%pearR
            #
            plt.axhline(linewidth=0.1)#, color="g")        # inc. width of x-axis and color it green
            plt.axvline(linewidth=0.1)#, color="r")        # inc. width of y-axis and color it red
            #disable tick labels
            frame1 = plt.gca()
            #frame1.axes.get_xaxis().set_visible(False)
            #frame1.axes.get_yaxis().set_visible(False)
            for xlabel_i in frame1.axes.get_xticklabels():
                #xlabel_i.set_visible(False)
                xlabel_i.set_fontsize(3)
            for xlabel_i in frame1.axes.get_yticklabels():
                xlabel_i.set_fontsize(3)
                #xlabel_i.set_visible(False)
            #logscale
            if logscale:
                plt.xscale('log')
                plt.yscale('log')
        #save title of subplot
        plt.title(title.strip(),fontsize=5)        
    #save title for entire fig
    pearR,pval = pearsonr(X,Y)
    plt.suptitle("%s-vs-%s r=%.5f" % (",".join(f.name.split(".")[0] for f in xbeds),",".join(f.name.split(".")[0] for f in ybeds),pearR))        
    #save or plot
    if outname:
        #create outdir if specified
        outdir = os.path.dirname(outname)
        if outdir and not os.path.isdir(outdir):
            os.makedirs(outdir)
        #save fig
        format = outname.split(".")[-1]
        plt.savefig(outname, dpi=dpi, facecolor='w', edgecolor='w',
          orientation='landscape', format=format, transparent=False)
    else:
        #plot fig
        plt.show()

def main():
    usage   = "usage: %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )

    parser.add_argument("-v", dest="verbose",  default=False, action="store_true")
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-x", dest="xbed",    nargs="+", type=file, 
                        help="counts in BED     [mandatory]" )
    parser.add_argument("-y", dest="ybed",    nargs="+", type=file, 
                        help="counts in BED     [mandatory]" )
    parser.add_argument("-o", dest="outname",   default="", #type=argparse.FileType('w'), 
                        help="output fname      [don't save, just show plot]" )    
    parser.add_argument("-l", dest="logscale", default=False, action="store_true", 
                        help="use log scale     [%(default)s]" )
    parser.add_argument("-n", dest="normalise", default=False, action="store_true", 
                        help="sample size normalisation [%(default)s]" )
    parser.add_argument("-c", dest="countCol", default=9, type=int, 
                        help="counts in column  [%(default)s]" )
    parser.add_argument("-d", dest="dpi",      default=300, type=int, 
                        help="figure dpi        [%(default)s]" )

    o = parser.parse_args()
    
    if o.verbose: 
        sys.stderr.write( "Options: %s\n" % str(o) )

    #plot
    bedcount2plot(o.xbed, o.ybed, o.outname, o.logscale, o.normalise, o.countCol, o.dpi, o.verbose)

if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
    