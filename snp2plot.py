#!/usr/bin/env python2
"""Process SNPs report(s) and plot density on chromosomes given some window size (kb).
It ignores contigs without SNPs and having size smaller than twice of the window.

It can handle vcf and bam2snps.py output files.

Author:
l.p.pryszcz+git@gmail.com

Barcelona, 21/03/2012
"""

import os, sys
from optparse import OptionParser
from datetime import datetime
from genome_annotation import load_gtf,load_gff,get_contig2size #bin/python_modules
import numpy as np
import matplotlib.pyplot as plt

def save_plot( outfn,window,positions,genes,contigSize,contig="",ax=None,width=0.5 ):
    """Save plotted figure.
    Return ax, usefull so next plots will have same scale.
    """
    # init figure object
    fig = plt.figure()
    if not ax:
        ax1 = fig.add_subplot(111)
    else:
        ax1 = fig.add_subplot(111,sharex=ax)
    plt.title( "Density of SNPs in %s [%s kb windows]" % ( contig,window ) )

    counts,coding = [],[]
    for pos in positions:
        # add empty bars
        i = int(pos/window)
        while i>=len(counts):
            counts.append( 0 )
            coding.append( 0 )
        # add all count
        counts[i] += 1
        # add coding only if overlap with gene
        if filter( lambda x: x[0]<=pos*1000<=x[1],genes ):
            coding[i] += 1
    
    xticks = [ i*window+(1.0*window-width*window)/2 for i in range(len(counts)) ]
    # add all snps as grey bars
    ax1.bar( xticks,counts,width=width*window,facecolor='black',alpha=0.3,label="non-coding" )
    # add coding as red
    ax1.bar( xticks,coding,width=width*window,facecolor='black',alpha=0.7,label="coding" )
 
    # ax1.legend( markerscale=(0.25,0.25),shadow=True )
    # remove top and right spines
    for loc, spine in ax1.spines.iteritems():
        if loc in ['left','bottom']:
            spine.set_position(('outward',10)) # outward by 10 points
        elif loc in ['right','top']:
            spine.set_color('none') # don't draw spine
        else:
            raise ValueError('unknown spine location: %s'%loc)

    # turn off ticks where there is no spine
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')

    # describe x and y axes
    ax1.set_xlabel( "Chromosome position [kb]" )
    ax1.set_ylabel( "Density" )
    # save plot
    ext = outfn.split('.')[-1]
    fig.savefig( outfn,format=ext,dpi=500 )
    return ax1


def vcf2snps( fn ):
    """Parse vcf file"""
    contig2snps = {}
    for l in open( fn ):
        l = l.strip()
        if not l or l.startswith("#"):
            continue
        contig,pos = l.split('\t')[:2]
        pos = int(pos)
        if contig not in contig2snps:
            contig2snps[contig] = []

        #add snp position
        contig2snps[contig].append( pos/1000.0 )
    return contig2snps

def file2snps( fn ):
    """Parse my own file"""
    #HE605202:236	15	T	1.0000	28	A	1.0000	intergenic
    contig2snps = {}
    for l in open( fn ):
        l = l.strip()
        if not l or l.startswith("#"):
            continue
        contig,pos = l.split('\t')[0].split(":")
        pos = int(pos)
        if contig not in contig2snps:
            contig2snps[contig] = []

        #add snp position
        contig2snps[contig].append( pos/1000.0 )
    return contig2snps    
    
def snps2plot( fn,window,contig2gene,contig2size,outbase,splitFn,ext,verbose ):
    """
    """
    if ".vcf" in fn:
        contig2snps = vcf2snps( fn )
    else:
        contig2snps = file2snps( fn )
    
    #plot
    bfn = os.path.basename( fn )
    if splitFn:
        bfn = bfn.split(".")[0]
    outdir = os.path.join( outbase,bfn )
    
    if not os.path.isdir( outdir ):
        os.makedirs( outdir )

    ax = None
    for contig in sorted( contig2snps.keys(),key=lambda x: contig2size[x],reverse=True ):
        #skip if less that 10 SNPs or chr shorter than 10kb
        if not contig2snps[contig] or contig2size[contig]<2*window*1000:
            continue
        #get contig size in kb
        contigSize = contig2size[contig]/1000.0
        if verbose:
            print contig,contigSize,len(contig2snps[contig])
        #generate outfn
        outfn = os.path.join( outdir,contig+"."+ext )
        ax = save_plot( outfn,window,contig2snps[contig],contig2gene[contig],contigSize,contig,ax )
    del ax
        
def main():
    
    usage = "usage: %prog [options] vcf1 [ vcf2 ... vcfN ]" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True

    parser.add_option("-g", dest="gtf",
                      help="genome annotation" )
    parser.add_option("-f", dest="fasta",
                      help="genome fasta" )
    parser.add_option("-o", dest="outbase", default="plots",
                      help="output directory [%default]" )
    parser.add_option("-s", dest="splitFn",  default=False, action="store_true", 
                      help="split fname (sheet name) by dot")
    parser.add_option("-w", dest="window",   default=10, type=int,
                      help="window size in kb [%default]")
    parser.add_option("-p", dest="ext",   default="png",
                      help="Supported: emf, eps, pdf, png, ps, raw, rgba, svg, svgz [%default]")
    parser.add_option("-v", dest="verbose",  default=True, action="store_false")
    
    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "%s\nFiles to process: %s\n" % ( str(o),", ".join( args ) ) )

    #check if any input file
    if not args:
        parser.error( "At least one input file has to be specified!" )

    #check if files exists
    for fn in args:
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" )
    
    #load genome - in fact need only contig sizes
    contig2size = get_contig2size( o.fasta )
    #load gtf
    if o.gtf.endswith(".gff"):
        gene2position, contig2gene = load_gff( o.gtf )
    else:
        gene2position, contig2gene = load_gtf( o.gtf )
    
    #process vcf
    for fn in args:
        print fn
        snps2plot( fn,o.window,contig2gene,contig2size,o.outbase,o.splitFn,o.ext,o.verbose )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
