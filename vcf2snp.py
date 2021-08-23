#!/usr/bin/env python2
"""
Parses vcf file (SNPs/INDELs) and check whether:
-coding
-nonsynonymous
-compare with reference (optionally)
In addition, filter for min depth of coverage and whether reads aligned to both strands.
"""

import os, sys
from optparse import OptionParser
from datetime import datetime
from genome_annotation import load_gtf,load_gff,genome2dict,coding_snp_info # bin/python_modules

def filtering( info,minDepth,minFreq,bothStrands=False ):
    """Return False if snp doesn't pass filtering.
    """
    # DP=7;AF1=1;CI95=0.5,1;DP4=0,0,4,1;MQ=26;FQ=-42
    if not "DP4=" in info:
        #sys.stderr.write( "Warning: parse_vcf: filtering: No DP4 in info: %s\n" % info )
        return True
    # get dp4
    dp4=info.split("DP4=")[1].split(";")[0]
    dp4list = [ int(x) for x in dp4.split(',') ] # rF,rR,aF,aR
    # depth
    if sum( dp4list )<minDepth:
        return
    # minFreq
    if 1 - sum( dp4list[:2] )*1.0/sum( dp4list[2:]) < minFreq:
        return
    # bothStrands
    if bothStrands and 0 in dp4list[2:]:
        return
    return True

def process_alts( ref,alts,contig,pos,contig2position,gene2position,contig2fasta,l ):
    """
    """
    outline = ""
    for alt in alts:
        genes = filter( lambda x: x[0]<pos+1<x[1], contig2position[contig] )
        if genes:
            for start,stop,feature,geneid in genes:
                # check effect of mutations
                if len(ref)!=len(alt): # indel
                    if feature == 'gene':
                        contig,CDSs,strand,function = gene2position[geneid] # get exons
                        cds = filter( lambda x: x[0]<pos+1<x[1],CDSs ) # check if overlaping with indel
                        if cds and len(ref)%3 != len(alt)%3:
                            outline += "%s\texonic\t%s\tframe-shift\t\t\t\t\t\t%s\n" % ( l,geneid,function )
                        else:
                            outline += "%s\tintronic\t%s\t\t\t\t\t\t%s\n" % ( l,geneid,function )
                    else:
                        outline += "%s\tintergenic\n" % l
                elif feature == 'gene':
                    contig,CDSs,strand,function = gene2position[geneid]
                    outline += "%s\t%s\t%s\n" % ( l,coding_snp_info( contig2fasta[contig],geneid,CDSs,strand,ref,alt,pos ),function )
                else:
                    outline += "%s\t%s\n" % ( l,feature )
        else:
            outline += "%s\tintergenic\n" % ( l, )
    
    return outline

def parse_vcf( fpath,outfn,contig2position,gene2position,contig2fasta,minDepth,minFreq,indels,bothStrains=True ):
    """
    """
    # select input
    if not fpath: # stdin if not infile
        handle = sys.stdin
    else:         # file
        handle = open(fpath)
    # select output
    if outfn:    # write to file, if specified
        out1 = open( outfn,'w')
    else:           # write to stdout
        out1 = sys.stdout

    # parse vcf
    snpsCount = indelsCount = 0
    for l in handle:
        l = l.strip()
        if l.startswith("#") or not l:
            if   l.startswith("##"):
                continue
            elif l.startswith("#CHROM"):
                l+="\tSNP type\tgene\tAA type\tAA position\tposition in codon\tref codon\tref AA\talt codon\talt AA"
            out1.write( l+"\n" )
            continue
        ## CHROM POS ID REF ALT QUAL FILTER INFO FORMAT bam
        contig,pos,id,ref,alt,qual,fitler,info,format = l[:-1].split("\t")[:9]
        # skip indels if not requested
        if info.startswith("INDEL") and not indels:
            continue
        # filter
        if not filtering( info,minDepth,minFreq,bothStrains ):
            continue

        # count
        if info.startswith("INDEL"):
            indelsCount += 1
        else:
            snpsCount   += 1

        # check if sysnonymous
        pos             = int(pos)
    
        # check if in gene
        if not contig in contig2position:
            sys.stderr.write("Warining: Contig %s is not present in GTF!\n" % contig )
            continue
        
        alts = alt.split(',') # G,T
        outline = process_alts( ref,alts,contig,pos,contig2position,gene2position,contig2fasta,l )
        out1.write( outline ) 

    sys.stderr.write( "SNPs:\t%s\nINDELs:\t%s\n" % ( snpsCount,indelsCount ) )

def main():
    
    usage = "usage: %prog [options]" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) # allow_interspersed_args=True

    parser.add_option("-g", dest="gtf",
                      help="genome annotation gtf/gff [requires -f]" )
    parser.add_option("-f", dest="fasta",
                      help="genome fasta [can be gzipped]" )
    parser.add_option("-i", dest="fpath",
                      help="input file [stdin]")
    parser.add_option("-o", dest="outfn",
                      help="output fname [stdout]")
    parser.add_option("-d", dest="minDepth", default=10,  type=int,
                      help="minimal depth [%default]")
    parser.add_option("-m", dest="minFreq",  default=0.8, type=float,
                      help="min frequency of alternative base [%default]")
    parser.add_option("-n", dest="indels",   default=True, action="store_false", 
                      help="ignore indels")
    parser.add_option("-b", dest="bothStrands", default=True, action="store_false", 
                      help="report events confirmed by single strand algs")
    parser.add_option("-v", dest="verbose",  default=True, action="store_false")
    
    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "%s\n" % ( str(o), ) )

    ctg2cds,id2gene,ctg2seq = {},{},{}
    if o.gtf: # if annotation
        # load genome
        if not o.fasta: # fasta has to be provided
            parser.errer( "Fasta file (-f) is requeired!" )
        elif not os.path.isfile( o.fasta ):
            parser.error( "No such file: %s" % o.fasta )
        ctg2seq        = genome2dict( o.fasta )

        # load genome annotation
        if not os.path.isfile( o.gtf ): # check if correct file
            parser.error( "No such file: %s" % o.gtf )
        # load gtf/gff
        if o.gtf.endswith(".gff"):
            id2gene,ctg2cds = load_gff( o.gtf )
        else:
            id2gene,ctg2cds = load_gtf( o.gtf )
        if o.verbose:
            sys.stderr.write( "Loaded annotation of %s CDS from %s\n" % ( len(id2gene),o.gtf ) )

    # parse pileup
    parse_vcf( o.fpath,o.outfn,ctg2cds,id2gene,ctg2seq,o.minDepth,o.minFreq,o.indels,o.bothStrands )
  
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
