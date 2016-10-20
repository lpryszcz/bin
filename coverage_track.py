#!/usr/bin/env python
desc="""Generate nice graph with coverage for multiple samples.

First generate cov.bed files is:
for f in bowtie2/*.bam; do echo `date` $f; bedtools coverage -counts -abam $f -b ../ensembl/Fusarium_oxysporum.FO2.15.dna.toplevel.fa.no_gaps.1kb.bed > $f.cov.bed; done

TODOs:
Normalise by GC?

"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 1/10/2012
"""

import os, pysam, sys
from datetime import datetime
from optparse import OptionParser
from math     import log
from Bio      import SeqIO
from Bio.SeqUtils  import GC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors
from Bio.Graphics  import GenomeDiagram
from genome_annotation import get_contig2coverage
from genome_annotation import load_sgd_gff,parse_gtf

#Categorical 12-step scheme, after ColorBrewer 11-step Paired Scheme from http://geography.uoregon.edu/datagraphics/color_scales.htm#Categorical Color Schemes
COLORCATS = (    
    [1.0, 0.5, 0.0],
    [0.2, 1.0, 0.0],
    [0.1, 0.7, 1.0],
    [0.4, 0.3, 1.0],
    [0.9, 0.1, 0.2],
    [1.0, 1.0, 0.2],
    [1.0, 0.75, 0.5],
    [0.7, 1.0, 0.55],
    [0.65, 0.93, 1.0],
    [0.8, 0.75, 1.0],
    [1.0, 0.6, 0.75],
    [1.0, 1.0, 0.6],
             )

def load_counts_beds( fnames,window,minreads,verbose ):
    """Load counts for multiple intervals from BED files.
    Return dict of counts for intervals and exp counts for
    given window size.
    Skip windows smaller than 0.25 * window size.
    """
    if verbose:
        sys.stderr.write( "Loading counts for %s BED files...\n" % len(fnames) )
    bedsdict  = {}
    expcounts,fnames2 = [],[]
    i = 0
    for fn in fnames:
        rcount = gsize = 0
        for l in open( fn ):
            if l.startswith("#") or not l:
                continue
            contig,s,e,c = l.split('\t')
            s,e,c = int(s),int(e),float(c)
            #skip windows smaller than 0.25 * window size
            if e-s < 0.25*window:
                continue
            #prepare list for contig
            if contig not in bedsdict:
                #bedsdict[contig] = [ [] for ii in range( len(fnames) ) ]
                bedsdict[contig] = [ [] ]
            while len( bedsdict[contig] ) < i+1:
                bedsdict[contig].append( [] )
            #add data
            bedsdict[contig][i].append( (s,e,c) )
            #update gsize and number of reads
            gsize  += e-s
            rcount += c
        #skip if not enough reads in sample
        if rcount < minreads * 10**6:
            if verbose:
                sys.stderr.write( " %s has not enough aligned reads (%s). Skipped\n" % ( fn,rcount) )
            continue
        #calculate exp count
        if verbose:
            sys.stderr.write( " %s: %.3f M reads\n" % (fn,rcount/10.0**6) )
        expcount = 1.0 * rcount * window / gsize
        expcounts.append( expcount )
        fnames2.append( fn )
        i += 1
    #sort
    for contig in bedsdict:
        for bedlist in bedsdict[contig]:
            bedlist.sort()
    return bedsdict,expcounts,fnames2

def seq2gcgraph( seq,bedlist ):
    """Return tuple of windows starts and GC """
    gcgraph = [ (0,0),] # ( s,GC( seq[s:e] ) ) for s,e,c in bedlist ]
    for s,e,c in bedlist:
        seqslice = str(seq[s:e])
        seqslice = seqslice.replace("N","").replace("n","")
        gcgraph.append( ( s+1,GC( seqslice ) ) )
    return gcgraph

def bam2graph( bam,contig,window ):
    """ """
    #first get mean coverage
    c2cs     = get_contig2coverage( bam )
    gsize    = sum( [ s for c,s in c2cs.itervalues() ] )
    rcount   = sum( [ c for c,s in c2cs.itervalues() ] )
    expcount = 1.0 * rcount * window / gsize
    #generate graph data
    gdata = []
    clen  = c2cs[contig][1]
    sam = pysam.Samfile( bam )
    for i in xrange( 0,clen,window ):
        #log2 = 0
        c = sam.count( reference=contig,start=i,end=i+window )
        if c:
            #last window is shorter!
            if i+window>clen:
                expcount = 1.0 * rcount * (clen-i) / gsize
            log2 = log( c/expcount,2 )
            gdata.append( (i+1,log2) )
    return gdata,clen

def _get_color(d):
    ci = round(255*d/8)
    if ci > 255:
        ci = 255
    return ci
    
def bed2SeqFeature( bedlist,window,expcount,cnvlimit ):
    """Generate feature data."""
    gdata = []
    for s,e,c in bedlist:
        #assume not reads and dels
        log2 = float('-inf')
        #get log2 only if any reads        
        if c:
            #get exp count
            expcountlocal = expcount
            if e-s < window:
                #if region shorter than window, normalize expcount accordingly
                expcountlocal = 1.0 * (e-s) / window * expcount
            log2 = log( c / expcountlocal,2 )
        #store dels and dups
        sf = SeqFeature(FeatureLocation(s,e))
        color = 0
        '''
        if   log2 <= -cnvlimit:
            color = (0,0,_get_color(-log2))
        elif log2 >=  cnvlimit:
            color = (_get_color(log2),0,0) #rgb
        '''
        if   log2 <= -3*cnvlimit:
            color = colors.darkblue
        elif log2 <= -2*cnvlimit:
            color = colors.blue
        elif log2 <= -1*cnvlimit:
            color = colors.lightblue
        elif log2 >=  3*cnvlimit:
            color = colors.darkred
        elif log2 >=  2*cnvlimit:
            color = colors.red
        elif log2 >=  1*cnvlimit:
            color = colors.lightsalmon #'''
        if color:
            gdata.append((sf,color))
    return gdata
    
def record2graph( fnames,beds,r,expcounts,minlog,window,verbose ):
    """ """
    #create diagram
    gdd  = GenomeDiagram.Diagram() #GDDiagram(gb)
    #add annotation
    gdt1 = gdd.new_track( 1,greytrack=1,name="%s: Genes & GC" % r.id,height=2.0,scale_smalltick_interval=5*10**4,scale_largetick_interval=25*10**4,scale_smallticks=0.15,scale_largetick_labels=1,scale_fontangle=0 )
    gdt1.greytrack_fontcolor = colors.black
    gdfs = gdt1.new_set("feature")
    for feature in r.features:
        if feature.type == "CDS":
            gdfs.add_feature( feature,colour=colors.grey )
    #add GC
    gdgs = gdt1.new_set("graph")
    #get gc graph
    gcgr = seq2gcgraph( r,beds[0] )
    gdgg = gdgs.new_graph( gcgr,"GC content",style="line",colour=colors.blue,center=50 )
    clen = gdgg.range()[1]
    basei  = 2
    height = 1.0
    #add coverage tracks for each bam file
    i = 0
    for bed,fn,expcount in zip( beds,fnames,expcounts ):
        gdata = bed2SeqFeature( bed,window,expcount,minlog )#; print i,gdata
        #adjust font
        fsize = 8-len(beds)/5
        if fsize < 1: fsize = 1
        gdt   = gdd.new_track( i+basei,greytrack=1,greytrack_fontsize=fsize,name=fn,height=height,scale_smalltick_interval=5*10**4,scale_smallticks=0.15,scale_largetick_labels=0,scale_fontangle=0 )
        gdt.greytrack_fontcolor = colors.black
        gdgs  = gdt.new_set("feature")
        #add CNVs
        for cnv,color in gdata:
            gdgs.add_feature( cnv,colour=color )
        i += 1
        
    #write 
    xl=xr  = 0.05
    width  = 841.8897637795275
    height = 595.275590551181
    '''if clen<10.0**6:
        xr = 1.0 - clen / 10.0**6 * 0.95
    else:
        width = clen * width / 10.0**6
    if len(fnames)>12:
        height = len(fnames)/12.0 * height'''
    #draw        
    gdd.draw( format="linear",pagesize=(width,height),xl=xl,xr=xr,orientation="landscape",tracklines=0,fragments=1,circular=0,track_size=0.75 ) # ,pagesize="A3"
    return gdd

def get_records( genome,gff,verbose ):
    """Return list of records from genome with
    annotations from annotation file.
    GFF is 1-base; inclusive so
    start-1, end the same
    """
    if verbose:
        sys.stderr.write( "Generating annotations for genome...\n" )
    #load annotation
    if gff.endswith(".gff"):
        contig2coding,trans2exon = load_sgd_gff( gff )
    else:
        contig2coding,trans2exon = parse_gtf( gff )
    
    records = []
    for r in SeqIO.parse( open(genome),"fasta" ):
        t = "CDS"
        if r.id in contig2coding:
            for s,e,name,stnd,score in contig2coding[r.id]:
                strand = 1
                if stnd == "-":
                    strand = -1
                r.features.append( SeqFeature(FeatureLocation(s-1,e),type=t,strand=strand,id=name ) )
        elif verbose:
            sys.stderr.write( " Warning: no annotation for %s\n" % r.id )
        #add to list
        records.append( r )
    return records
    
def coverage_graph( outdir,fnames,gfn,fnformat,gff,lenlimit,ext,minlog,minreads,window,verbose ):
    """Generate plot for depth of coverage"""
    #load beds
    countsdict,expcounts,fnames = load_counts_beds( fnames,window,minreads,verbose )
    if not fnames:
        sys.stderr.write( "No files left after filtering!\n" )
        return
    elif len(fnames) == 1:
        sys.stderr.write( "One file left after filtering!\n" )
        return
    #create outdir
    if not os.path.isdir( outdir ):
        os.makedirs( outdir )
    if verbose:
        sys.stderr.write( "Generating genome graphs...\n" )
    #process chromosomes/contigs independently
    if fnformat in set( ("genbank","embl","gb") ):
        seqobjects = SeqIO.parse( open(gfn),fnformat )
    else:
        seqobjects = get_records( gfn,gff,verbose )
    for r in seqobjects:
        contig = r.name
        if len(r.seq) >= lenlimit*10**3:
            if verbose:
                sys.stderr.write( " %s          \r" % contig )
            #get graph
            counts = countsdict[contig]
            gdd   = record2graph( fnames,counts,r,expcounts,minlog,window,verbose )
            #save
            outfn = os.path.join( outdir,"%s.%s" % ( contig,ext ) )
            gdd.write( outfn,ext )

def main():
    usage   = "usage: %prog [options]"
    version = "%prog 1.0"
    parser  = OptionParser( usage=usage,version=version,description=desc,epilog=epilog ) #allow_interspersed_args=True

    parser.add_option("-v", dest="verbose",  default=False, action="store_true")
    parser.add_option("-i", dest="genome",   default="", 
                      help="genome file         [%default]" )
    parser.add_option("-f", dest="format",   default="genbank", 
                      help="genome file format  [%default]\nAllowed: genbank, embl" )    
    parser.add_option("-g", dest="gff",      default="", 
                      help="annotation file     [required if no gb/embl genome]" )
    parser.add_option("-o", dest="outdir",   default="coverage_plot", 
                      help="output directory    [%default]" )
    parser.add_option("-e", dest="ext",      default="pdf", 
                      help="outfile extension   [%default]" )
    parser.add_option("-l", dest="lenlimit", default=100, type=int, 
                      help="only chr above l kb [%default kb]" )
    parser.add_option("-m", dest="mincnv",   default=0.5, type=float,
                      help="log2 to call CNV    [%default]" )
    parser.add_option("-r", dest="minreads", default=5, type=int,
                      help="min r million reads [%default millions]" )
    parser.add_option("-w", dest="window",   default=1000, type=int, 
                      help="window size         [%default bp]" )
    
    o,args = parser.parse_args()
    if o.verbose: 
        sys.stderr.write( "Options: %s\nArgs: %s\n" % (str(o),", ".join(args)) )

    if not args:
        parser.error( "Provide at least one BED file!" )
    
    for fn in args:
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )

    if not o.genome:
        parser.error( "Genome file has to be specified" )
    if not os.path.isfile( o.genome ):
        parser.error( "No such file: %s" % o.genome )

    if o.format not in set(["genbank","embl","gb"]) and not o.gff:
        parser.error( "Specify annotation file (gff)!" ) #"Only genbank/embl genome files are accepted" )
      
    #plot
    coverage_graph( o.outdir,args,o.genome,o.format,o.gff,o.lenlimit,o.ext,o.mincnv,o.minreads,o.window,o.verbose )
  
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )