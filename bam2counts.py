#!/usr/bin/env python
"""Return number of reads aligned to each gene.

Author:
l.p.pryszcz@gmail.com

Barcelona, 25/04/2012
"""

import commands, math, os, sys
from optparse import OptionParser
from datetime import datetime
from genome_annotation import load_gtf, load_gff, genome2dict, get_gc

def _count( bam,contig,start,end ):
    """Return number of reads from given region"""
    cmd   = "samtools view -c %s %s:%s-%s" % ( bam,contig,start,end )
    count = commands.getoutput( cmd )
    if not count.isdigit():
        sys.stderr.write( " Error: samtools: %s\n" % count )
        sys.exit()
    count = int( count )
    return count
    

def bam2counts( bam,normalise,id2gene,ctg2cds,ctg2seq,verbose ):
    """Print counts for every gene.
    RPKM normalisation if resequsted."""
    # get readcount
    if verbose:
        sys.stderr.write( "Counting reads on genes...\n" )
    i = totReads = 0
    gene2counts = {}
    for contig in ctg2cds:#.keys()[:1]:
        i+=1
        if verbose:
            sys.stderr.write( " %s ( %s / %s )          \r" % ( contig,i,len(ctg2cds) ) )

        for start,end,feature,geneid in ctg2cds[contig]:
            count     = _count( bam,contig,start,end )
            totReads += count

            strand    = id2gene[geneid][2]
            gene2counts[geneid] = ( contig,start,end,count )
            
    #print info
    if normalise:
        if verbose:
            sys.stderr.write( "Counting reads on exons...\n" )
        sys.stdout.write( 'chro\tstrand\tstart\tend\ttype\ttranscript\tgene_symbol\tlength\tcoverage\tRPKM\tlog2(RPKM+1)\tgc\n' )
    for geneid in sorted( gene2counts.keys() ):
        contig, start, end, count = gene2counts[geneid]
        line = "%s\t%s\n" % (geneid,count)
        if normalise: #reads per kb of gene per million of aligned reads
            glen   = end - start
            rpkm   = count * 10.0**9 / totReads / glen
            gc     = get_gc( ctg2seq,contig,start,end )#; print contig,start,end,gc            
            strand = id2gene[geneid][2]
            #chr strand start end transcript/exon exon_id transcript_id length coverage rpkm log2(rpkm+1) gc
            line   = "%s\t%s\t%s\t%s\ttranscript\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ( contig,strand,start,end,geneid,geneid,glen,float(count)/glen,rpkm,math.log(rpkm+1,2),gc )

            #add exons
            ti = 0
            for start,end in id2gene[geneid][1]:
                ti   += 1
                count = _count( bam,contig,start,end )
                glen  = end - start
                if not glen:
                    sys.stderr.write(" Error: Zero exon length @ %s %s:%s-%s\n" % (geneid,contig,start,end) )
                    continue
                rpkm  = count * 10.0**9 / totReads / glen
                gc    = get_gc( ctg2seq,contig,start,end )#; print contig,start,end,gc            

                line += "%s\t%s\t%s\t%s\texon\t%s.%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ( contig,strand,start,end,geneid,ti,geneid,glen,float(count)/glen,rpkm,math.log(rpkm+1,2),gc )

        sys.stdout.write( line )
    
def main():

    usage  = "usage: %prog [options]\nfor f in *.bam; do echo `date` $f; bam2counts.py -rv -i $f -g F.oxysporum.gtf > $f.genecounts.txt; done" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True

    parser.add_option("-i", dest="bam", default="",
                      help="bam file")
    parser.add_option("-g", dest="gtf",default="",
                      help="genome annotation gtf/gff" )
    parser.add_option("-r", dest="rpkm",  default=False, action="store_true",
                      help="RPKM normalisation (reads per kb of gene per million of aligned reads)" )
    parser.add_option("-f", dest="fasta", default="",
                      help="genome fasta [required if -r]")
    parser.add_option("-v", dest="verbose",  default=False, action="store_true" )    
  
    ( o, fnames ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nArgs: %s\n" % ( o,fnames ) )

    for fn in ( o.bam,o.gtf ):
        if not fn:
            parser.error( "Provide input file!" )
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )
        
    ctg2cds,id2gene,ctg2seq = {},{},{}
    # load gtf/gff
    if o.gtf:
        if o.gtf.endswith(".gff"):      
            id2gene,ctg2cds = load_gff( o.gtf )
        else:
            id2gene,ctg2cds = load_gtf( o.gtf )        
    if o.verbose:
        sys.stderr.write( "Loaded annotation of %s CDS from %s\n" % ( len(id2gene),o.gtf ) )

    if o.rpkm:
        if not o.fasta:
            parser.error( "Specify genome fasta file!" )
        if not os.path.isfile( o.fasta ):
            parser.error( "No such file: %s" % o.fasta )
        ctg2seq = genome2dict( o.fasta )

    bam2counts( o.bam,o.rpkm,id2gene,ctg2cds,ctg2seq,o.verbose )
    
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
