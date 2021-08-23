#!/usr/bin/env python2
desc="""Annotate SNPs with gene(s) and functions.

CHANGELOG:
+ 1.1:
- heterozygous sites support

"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 28/06/2012
"""

import argparse, os, sys
from Bio      import SeqIO
from datetime import datetime
from genome_annotation import load_gtf, load_gff, genome2dict, coding_snp_info 

def load_tab(tabfn):
    """Return gene/transcript annotation from TAB file ie interpor"""
    gene2fun = {}
    for l in open(tabfn):
        gene, function = l[:-1].split('\t')[:2]
        if gene not in gene2fun:
            gene2fun[gene] = []
        gene2fun[gene].append(function)
    return gene2fun

def load_pfam(pfam):
    """Return gene/transcript to pfam matches dictionary."""
    sys.path.append("/users/tg/lpryszcz/bin")
    import pfam2gtf
    pfam_data = pfam2gtf.load_gene2pfam(pfam)
    gene2pfam = {}
    for gene in pfam_data:
        pfams = ""
        for name,acc,desc in pfam_data[gene]:
            pfams += "%s [%s]; " % ( desc,acc )
        gene2pfam[gene] = pfams[:-2]
    return gene2pfam

def load_fasta_headers(fasta):
    """Return gene/transcript to header."""
    gene2ann = {}
    for r in SeqIO.parse( open(fasta),"fasta" ):
        # >FOXG_00006T0 | FOXG_00006 | Fusarium oxysporum f. sp. lycopersici 4287 protein phosphatase PP2A regulatory subunit B (475 aa)
        gene = r.id
        gene2ann[gene] = " ".join( r.description.split("|")[-1].split()[6:] )

    return gene2ann

def process_alt(ref, alt, contig, pos, contig2position, gene2position, contig2fasta, \
                trans2ann, trans2pfam, trans2tab, l):
    """ """
    outline = ""
    genes = []
    if contig in contig2position:
        genes = filter(lambda x: x[0]<pos+1<x[1], contig2position[contig])
    if genes:
        for start,stop,feature,geneid in genes:
            if feature == 'gene':
                fastaAnn = pfamAnn = tabAnn = ""
                if geneid in trans2ann:
                    fastaAnn = trans2ann[geneid]
                if geneid in trans2pfam:
                    pfamAnn  = trans2pfam[geneid]
                if geneid in trans2tab:
                    tabAnn   = trans2tab[geneid]
                contig,CDSs,strand,function,frame = gene2position[geneid]
                outline += "%s\t%s\t%s\t%s\t%s\t%s\n" % (l, coding_snp_info(contig2fasta[contig], geneid, CDSs, strand, ref, alt, pos), function, fastaAnn, pfamAnn, "; ".join(tabAnn))
            else:
                outline += "%s\t%s\n" % (l, feature)
    else:
        outline += "%s\tintergenic\n" % (l,)
    
    return outline

def parse_snps(handle, out, contig2position, gene2position, contig2fasta, \
               trans2ann, trans2pfam, trans2tab, verbose):
    """Process SNPs."""
    # change reference
    if verbose:
        sys.stderr.write("Adjusting reference sequence...\n")
    for l in handle:
        l = l.strip()
        if l.startswith("#") or not l:
            continue
        ##
        lData = l[:-1].split("\t")
        if len(lData)<5:
            if verbose:
                sys.stderr.write("[WARNING] Wrong line: %s\n"%str(lData))
            continue
        if ':' in lData[0]:
            #coord,refCov,refBase,altBase = lData[:4]
            coord,refCov,refBase,refFreq,altCov,altBase,altFreq = lData[:7]
            #alter contig2fasta
            contig,pos = coord.split(':')
        else:
            contig,pos,vid,refBase,altBase = lData[:5]
        pos        = int(pos)
        #unload alternatives 
        refBases = set(refBase.split(','))
        altBases = set(altBase.split(','))
        if len(refBases)>len(altBases):
            dbases = refBases.difference(altBases)
            if len(dbases)>1:
                sys.stderr.write("[Warning] Including 1 of %s reference bases\n"%len(dbases))
            #print refBases, altBase, dbases
            refBase = dbases.pop()
        else:
            refBase = refBases.pop()
            
        if contig not in contig2fasta:
            if verbose:
                sys.exit("Warning: %s not in genome\n" % contig )

        if contig2fasta[contig][pos-1] != refBase:
            if verbose:
                sys.stderr.write( " %s %s > %s\n" % ( coord,contig2fasta[contig][pos],refBase ) )
            contig2fasta[contig] = contig2fasta[contig][:pos-1] + refBase + contig2fasta[contig][pos:]
    
    # parse snp file
    snpsCount = indelsCount = 0
    headeradded = 0
    pcontig = ppos = 0
    storage = []
    handle.seek(0)
    for l in handle: #open(fpath):
        l = l.strip()
        if l.startswith("#") or not l:
            if   l.startswith("##") or headeradded:
                continue
            elif l.startswith("#"):
                headeradded = 1 
                l+="\tSNP type\tgene\tAA type\tAA position\tposition in codon\tref codon\tref AA\talt codon\talt AA\tfuntcion\tfasta annotation\tpfam\ttab annotation\n"
            out.write(l)
            continue
        ##
        snpsCount += 1
        lData = l[:-1].split("\t")
        if len(lData)<5:
            if verbose:
                sys.stderr.write("[WARNING] Wrong line: %s\n"%str(lData))
            continue
        if ':' in lData[0]:
            coord, refCov, refBase, refFreq, altCov, altBase, altFreq = lData[:7]
            #alter contig2fasta
            contig, pos = coord.split(':')
        else:
            contig, pos, vid, refBase, altBase = lData[:5]
        pos        = int(pos)
        #unload alternatives
        refBases = set(refBase.split(','))
        altBases = set(altBase.split(','))
        added = set()
        for refBase in refBases.difference(altBases):
            for altBase in filter(lambda x: x != refBase, altBases): 
                outline = process_alt(refBase, altBase, contig, pos, contig2position, gene2position, contig2fasta, trans2ann, trans2pfam, trans2tab, l)
                out.write(outline)
                added.add((refBase,altBase))
        for altBase in altBases.difference(refBases):
            for refBase in filter(lambda x: x != altBase, refBases):
                if (refBase,altBase) in added:
                    continue
                outline = process_alt(refBase, altBase, contig, pos, contig2position, gene2position, contig2fasta, trans2ann, trans2pfam, trans2tab, l)
                out.write(outline)
        
 
    #sys.stderr.write( "SNPs:\t%s\nINDELs:\t%s\n" % ( snpsCount,indelsCount ) )
    sys.stderr.write("SNPs:\t%s\n"%(snpsCount,))

def main():
    
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.1')
    parser.add_argument("-g", "--gtf",
                        help="genome annotation gtf/gff [requires -f]" )
    parser.add_argument("-f", "--fasta",
                        help="genome fasta [can be gzipped]" )
    parser.add_argument("-i", "--input", type=file,  #default=sys.stdin,
                        help="input stream [stdin]")
    parser.add_argument("-o", "--out", default=sys.stdout, 
                        help="output stream [stdout]")
    parser.add_argument("-p", "--pfam", default="", 
                        help="pfam tblout file") 
    parser.add_argument("-q", "--faa", default="", 
                        help="proteome fasta (to get protein annotation)") 
    parser.add_argument("-t", "--tab", default="", 
                        help="tab-delimited annotation") 
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    ctg2cds, id2gene, ctg2seq = {},{},{}
    if o.gtf: # if annotation
        # load genome
        if not o.fasta: # fasta has to be provided
            parser.errer("Fasta file (-f) is requeired!")
        elif not os.path.isfile( o.fasta ):
            parser.error("No such file: %s"%o.fasta)
        ctg2seq        = genome2dict(o.fasta)

        # load genome annotation
        if not os.path.isfile(o.gtf): # check if correct file
            parser.error("No such file: %s"%o.gtf)
        # load gtf/gff
        if o.gtf.endswith(".gff"):
            id2gene,ctg2cds = load_gff(o.gtf)
        else:
            id2gene,ctg2cds = load_gtf(o.gtf)
        if o.verbose:
            sys.stderr.write("Loaded annotation of %s CDS from %s\n"%(len(id2gene), o.gtf))

    #load function annotation
    trans2ann = trans2pfam = trans2tab = {}
    if o.faa:
        trans2ann = load_fasta_headers(o.faa)
    if o.pfam:
        trans2pfam = load_pfam(o.pfam)
    if o.tab:
        trans2tab = load_tab(o.tab)
    # parse pileup
    parse_snps(o.input, o.out, ctg2cds, id2gene, ctg2seq, trans2ann, trans2pfam, \
               trans2tab, o.verbose)
  
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
