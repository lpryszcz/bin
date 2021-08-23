#!/usr/bin/env python2
desc="""Generate nice graph with coverage for multiple samples.
Tracks are in colors and merged track is added if less than 13 samples. 

First generate cov.bed files is:
for f in bowtie2/*.bam; do echo `date` $f; bedtools coverage -counts -abam $f -b ../ensembl/Fusarium_oxysporum.FO2.15.dna.toplevel.fa.no_gaps.1kb.bed > $f.cov.bed; done

TODOs:
Normalise by GC?

"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 1/10/2012
"""

import argparse, os, pysam, sys
from datetime import datetime
#from optparse import OptionParser
from math     import log
from Bio      import SeqIO
from Bio.SeqUtils  import GC
from reportlab.lib import colors
from Bio.Graphics  import GenomeDiagram
from genome_annotation import get_contig2coverage,load_sgd_gff,parse_gtf
from Bio.SeqFeature import SeqFeature, FeatureLocation

def load_counts_beds( files,window,minreads,verbose ):
    """Load counts for multiple intervals from BED files.
    Return dict of counts for intervals and exp counts for
    given window size.
    Skip windows smaller than 0.25 * window size.
    """
    if verbose:
        sys.stderr.write( "Loading counts for %s BED files...\n" % len(files) )
    bedsdict  = {}
    expcounts,fnames = [],[]
    i = 0
    for f in files:
        rcount = gsize = 0
        for l in f:
            if l.startswith("#") or not l:
                continue
            contig,s,e,c = l.split('\t')
            s,e,c = int(s),int(e),int(c)
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
        '''#skip if not enough reads in sample
        if rcount < minreads * 10**6:
            if verbose:
                sys.stderr.write( " %s has not enough aligned reads (%s). Skipped\n" % ( fn,rcount) )
            continue'''
        #calculate exp count
        if verbose:
            sys.stderr.write( " %s: %s entries\n" % (f.name,'{0:,}'.format(rcount) ) )
        expcount = 1.0 * rcount * window / gsize
        expcounts.append( expcount )
        fnames.append( f.name )
        i += 1
    #sort
    for contig in bedsdict:
        for bedlist in bedsdict[contig]:
            bedlist.sort()
    return bedsdict,expcounts,fnames

def seq2gcgraph( seq,bedlist ):
    """Return tuple of windows starts and GC """
    gcgraph = [ (0,0),] # ( s,GC( seq[s:e] ) ) for s,e,c in bedlist ]
    for s,e,c in bedlist:
        seqslice = str(seq[s:e])
        seqslice = seqslice.replace("N","").replace("n","")
        gcgraph.append( ( s+1,GC( seqslice ) ) )
    return gcgraph

def bed2graph( bedlist,window,expcount,minlog ):
    """Generate graph data."""
    gdata = [ (0,minlog) ]
    for s,e,c in bedlist:
        log2 = get_log2(s,e,c,expcount,window,minlog)
        #don't allow too small values of log2 
        if log2<minlog:
            log2 = minlog
        if log2>-minlog:
            log2 = -minlog
        gdata.append( (s+1,log2) )
    return gdata

def get_log2(s,e,c,expc,window,minlog=float('-inf')):
    """ """
    #assume not reads and dels
    log2 = minlog
    #get log2 only if any reads        
    if c:
        #get exp count
        expcountlocal = expc
        if e-s < window:
            #if region shorter than window, normalize expcount accordingly
            expcountlocal = 1.0 * (e-s) / window * expc
        log2 = log( c / expcountlocal,2 )
    return log2
    
def bed2SeqFeature(bed1, expcount1, bed2, expcount2, window):
    """Generate feature data:
    - mark heterozygous SNP-rich with orange
    - mark homozygous SNP-rich with darkgrey
    """
    gdata = []
    bed1.sort()
    bed2.sort()
    for (s,e,c1),(s2,e2,c2) in zip(bed1,bed2):
        # hapA
        color = "darkgrey"#0
        #mark heterozygous SNP-rich with orange # hetero
        if   get_log2(s,e,c2,expcount2,window)>0:
            color = 0 #"orange"
        #mark homozygous SNP-rich with darkgrey # hapB
        elif get_log2(s,e,c1,expcount1,window)>0:
            color = "orange" #"darkgrey"
            
        #store dels and dups
        if color:
            sf = SeqFeature(FeatureLocation(s,e))
            gdata.append((sf,color))
    return gdata
    
def record2graph( fnames,beds,r,minlog,window,verbose ):
    """ """
    #create diagram
    gdd  = GenomeDiagram.Diagram() #GDDiagram(gb)
    #add annotation
    gdt1 = gdd.new_track( 1,greytrack=1,name="Genes and GC",height=1.5,
                          scale_smalltick_interval=5*10**4,scale_smallticks=0.15,
                          scale_largeticks=1.0,scale_largetick_labels=1,scale_largetick_interval=250*10**3,
                          scale_fontangle=0 )
    gdt1.greytrack_fontcolor = colors.black
    gdfs = gdt1.new_set("feature")
    for feature in r.features:
        if feature.type == "CDS":
            gdfs.add_feature( feature,color=colors.grey )
    #add GC
    gdgs = gdt1.new_set("graph")
    #get gc graph
    gcgraph = seq2gcgraph( r,beds[0][0][0] )
    gdgs.new_graph( gcgraph,"GC content",style="line",color=colors.blue,center=50 )

    #add coverage tracks for each bed tracks file
    ggraphs,gtracks = [],[]
    for i,fn in enumerate( fnames ):
        #add individual track
        gdt   = gdd.new_track( i+2,greytrack=1,name=fn,height=1,scale_smalltick_interval=5*10**4,scale_smallticks=0.15,scale_largetick_labels=0,scale_fontangle=0 )
        gdt.greytrack_fontcolor = colors.black
        #adjust font
        fsize = 8-len(fnames)/5
        if fsize < 1:
            fsize = 1        
        gdt.greytrack_fontsize = fsize
        #add feature and graph
        gtracks.append( gdt.new_set("feature") )
        ggraphs.append( gdt.new_set("graph") )
        
    j=0
    #colorstuple = ("blue","darkgrey","orange")
    #first add coverage track
    i=0
    for bed,expcount,bed1,expcount1,bed2,expcount2 in zip(beds[0][0], beds[0][1], beds[1][0], beds[1][1], beds[2][0], beds[2][1]):
        #add SNPs colored data
        fdata = bed2SeqFeature( bed1,expcount1,bed2,expcount2,window )
        for feature,color in fdata:
            gtracks[i].add_feature( feature,colour=color )
        #add coverage track
        gdata = bed2graph( bed,window,expcount,minlog )
        gdgg  = ggraphs[i].new_graph( gdata,fn,style="line",linewidth=0.5,center=0,color="blue" )
            
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
    
def coverage_graph( outdir,covbedfiles,homos,heteros,genome,fnformat,gff,lenlimit,ext,minlog,minreads,window,verbose ):
    """Generate plot for depth of coverage"""
    #load beds
    countsdict,expcounts,fnames = load_counts_beds( covbedfiles,window,minreads,verbose )
    if not fnames:
        sys.stderr.write( "No files left after filtering!\n" )
        return
    elif len(fnames) == 1:
        sys.stderr.write( "One file left after filtering!\n" )
        return
    #load homo and heterosnps
    if homos:
        homocountsdict,homoexpcounts,homofns = load_counts_beds( homos,window,0,verbose )
    if heteros:
        hetecountsdict,heteexpcounts,hetefns = load_counts_beds( heteros,window,0,verbose )
    #create outdir
    if not os.path.isdir( outdir ):
        os.makedirs( outdir )
    if verbose:
        sys.stderr.write( "Generating genome graphs...\n" )
    #process chromosomes/contigs independently
    if fnformat in set( ("genbank","embl","gb") ):
        seqobjects = SeqIO.parse( genome,fnformat )
    else:
        seqobjects = get_records( genome.name,gff,verbose )        
    for r in seqobjects:
        contig = r.name
        if len(r.seq) >= lenlimit*10**3:
            if verbose:
                sys.stderr.write( " %s          \r" % contig )
            
            if homos:
                homostuple = (homocountsdict[contig],homoexpcounts)
            else:
                homostuple = ([],[])
            if heteros:
                heterostuple = (hetecountsdict[contig],heteexpcounts)
            else:
                heterostuple = ([],[])
            #get graph    
            counts = [ (countsdict[contig],expcounts),homostuple,heterostuple ]
            gdd    = record2graph( fnames,counts,r,minlog,window,verbose )
            #save
            outfn  = os.path.join( outdir,"%s.%s" % ( contig,ext ) )
            gdd.write( outfn,ext )

def main():
    usage   = "usage: %(prog)s [options]"
    #version = "%prog 1.0"
    #parser  = OptionParser( usage=usage,version=version,description=desc,epilog=epilog ) #allow_interspersed_args=True
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )

    parser.add_argument("-v", dest="verbose",  default=False, action="store_true")
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="genome",   type=file, #argparse.FileType('r'),  
                        help="genome file         [%(default)s]" )
    parser.add_argument("-f", dest="format",   default="genbank", #type=complex, choices=["genbank","embl"], 
                        help="genome file format  [%(default)s] genbank or embl" )    
    parser.add_argument("-g", dest="gff",      default="", 
                        help="annotation file     [required if no gb/embl genome]" )
    parser.add_argument("-o", dest="outdir",   default="genome_plot", 
                        help="output directory    [%(default)s]" )
    parser.add_argument("-c", dest="covbed",   nargs='+', type=file,
                        help="coverage BED files  [%(default)s]" )
    parser.add_argument("-s", dest="homosnp",  nargs='+', type=file,
                        help="Homozygous SNP BED files [%(default)s]" )
    parser.add_argument("-t", dest="hetesnp",  nargs='+', type=file,
                        help="Heterozygous SNP BEDs files [%(default)s]" )
    parser.add_argument("-e", dest="ext",      default="pdf", 
                        help="outfile extension   [%(default)s]" )
    parser.add_argument("-l", dest="lenlimit", default=10, type=int, 
                        help="only chr above l kb [%(default)s kb]" )
    parser.add_argument("-m", dest="minlog",   default=-4, type=int,
                        help="min log2 value      [%(default)s kb]" )
    parser.add_argument("-r", dest="minreads", default=5, type=int,
                        help="min r million reads [%(default)s millions]" )
    parser.add_argument("-w", dest="window",   default=1000, type=int, 
                        help="window size         [%(default)s bp]" )
    
    o = parser.parse_args()
    
    if o.verbose: 
        sys.stderr.write( "Options: %s\n" % str(o) )

    if o.format not in set(["genbank","embl","gb"]) and not o.gff:
        parser.error( "Specify annotation file (gff)!" ) #"Only genbank/embl genome files are accepted" )
    #plot
    coverage_graph( o.outdir,o.covbed,o.homosnp,o.hetesnp,o.genome,o.format,o.gff,o.lenlimit,o.ext,o.minlog,o.minreads,o.window,o.verbose )
  
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )