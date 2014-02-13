#!/usr/bin/env python
"""Annotate SNPs with gene(s) and functions
"""

import os, sys
from Bio      import SeqIO
from optparse import OptionParser
from datetime import datetime
from genome_annotation import load_gtf,load_gff,genome2dict,coding_snp_info # bin/python_modules

def load_pfam( pfam ):
    """Return gene/transcript to pfam matches dictionary."""
    sys.path.append( "/users/tg/lpryszcz/bin" )
    import pfam2gtf
    pfam_data = pfam2gtf.load_gene2pfam( pfam )
    gene2pfam = {}
    for gene in pfam_data:
        pfams = ""
        for name,acc,desc in pfam_data[gene]:
            pfams += "%s [%s]; " % ( desc,acc )
        gene2pfam[gene] = pfams[:-2]

    return gene2pfam

def load_fasta_headers( fasta ):
    """Return gene/transcript to header."""
    gene2ann = {}
    for r in SeqIO.parse( open(fasta),"fasta" ):
        # >FOXG_00006T0 | FOXG_00006 | Fusarium oxysporum f. sp. lycopersici 4287 protein phosphatase PP2A regulatory subunit B (475 aa)
        gene = r.id
        gene2ann[gene] = " ".join( r.description.split("|")[-1].split()[6:] )

    return gene2ann

def process_alt( ref,alt,contig,pos,contig2position,gene2position,contig2fasta,trans2ann,trans2pfam,l ):
    """
    """
    outline = ""

    genes = filter( lambda x: x[0]<pos+1<x[1], contig2position[contig] )
    if genes:
        for start,stop,feature,geneid in genes:
            if feature == 'gene':
                fastaAnn = pfamAnn = ""
                if geneid in trans2ann:
                    fastaAnn = trans2ann[geneid]
                if geneid in trans2pfam:
                    pfamAnn  = trans2pfam[geneid]
                contig,CDSs,strand,function,frame = gene2position[geneid]
                outline += "%s\t%s\t%s\t%s\t%s\n" % (l, coding_snp_info(contig2fasta[contig], geneid, CDSs, strand, ref, alt, pos), function, fastaAnn, pfamAnn)
            else:
                outline += "%s\t%s\n" % (l, feature)
    else:
        outline += "%s\tintergenic\n" % (l,)
    
    return outline

def parse_snps( fpath,outfn,contig2position,gene2position,contig2fasta,trans2ann,trans2pfam,verbose ):
    """
    """
    # select output
    if outfn:    # write to file, if specified
        out1 = open( outfn,'w')
    else:           # write to stdout
        out1 = sys.stdout

    # change reference
    if verbose:
        sys.stderr.write( "Adjusting reference sequence...\n" )
    for l in open(fpath):
        l = l.strip()
        if l.startswith("#") or not l:
            continue
        ##
        lData = l[:-1].split("\t")
        if ':' in lData[0]:
            #coord,refCov,refBase,altBase = lData[:4]
            coord,refCov,refBase,refFreq,altCov,altBase,altFreq = lData[:7]
            #alter contig2fasta
            contig,pos = coord.split(':')
        else:
            contig,pos,vid,refBase,altBase = lData[:5]
        pos        = int(pos)

        if contig not in contig2fasta:
            sys.exit("Warning: %s not in genome\n" % contig )

        if contig2fasta[contig][pos-1] != refBase:
            if verbose: sys.stderr.write( " %s %s > %s\n" % ( coord,contig2fasta[contig][pos],refBase ) )
            contig2fasta[contig] = contig2fasta[contig][:pos-1] + refBase + contig2fasta[contig][pos:]
        
    # parse snp file
    snpsCount = indelsCount = 0
    headeradded = 0
    for l in open(fpath):
        l = l.strip()
        if l.startswith("#") or not l:
            if   l.startswith("##") or headeradded:
                continue
            elif l.startswith("#"):
                headeradded = 1 
                l+="\tSNP type\tgene\tAA type\tAA position\tposition in codon\tref codon\tref AA\talt codon\talt AA\tfuntcion\tfasta annotation\tpfam\n"
            out1.write( l )
            continue
        ##
        snpsCount += 1
        lData = l[:-1].split("\t")
        if ':' in lData[0]:
            coord,refCov,refBase,refFreq,altCov,altBase,altFreq = lData[:7]
            #alter contig2fasta
            contig,pos = coord.split(':')
        else:
            contig,pos,vid,refBase,altBase = lData[:5]
        pos        = int(pos)
        
        # check if in gene
        if not contig in contig2position:
            sys.stderr.write("Warining: Contig %s is not present in GTF!\n" % contig )
            continue
        
        outline = process_alt( refBase,altBase,contig,pos,contig2position,gene2position,contig2fasta,trans2ann,trans2pfam,l )
        out1.write( outline ) 

    #sys.stderr.write( "SNPs:\t%s\nINDELs:\t%s\n" % ( snpsCount,indelsCount ) )
    sys.stderr.write( "SNPs:\t%s\n" % ( snpsCount, ) )

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
    parser.add_option("-p", dest="pfam", default="", 
                      help="pfam tblout file") 
    parser.add_option("-q", dest="faa", default="", 
                      help="proteome fasta (to get protein annotation)") 
    parser.add_option("-v", dest="verbose",  default=False, action="store_true")
    
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

    #load function annotation
    trans2ann = trans2pfam = {}
    if o.faa:
        trans2ann = load_fasta_headers( o.faa )
    if o.pfam:
        trans2pfam = load_pfam( o.pfam )
            
    # parse pileup
    parse_snps( o.fpath,o.outfn,ctg2cds,id2gene,ctg2seq,trans2ann,trans2pfam,o.verbose )
  
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
