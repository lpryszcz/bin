#!/usr/bin/env python
desc="""Parse augustus GFF and save genes, cds and protein sequences.
Not reporting flanking sequence so far.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 1/08/2013
"""

import argparse, os, sys
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from genome_annotation import load_gtf, genome2dict

def get_cds_pep(contigSeq, gene, boundaries, strand, function, frames, codonTable):
    """Return CDS and amino acid seq"""
    #concatenate CDS seq
    seq = Seq("")
    for s, e in boundaries:
        seq += contigSeq[s-1:e]
    #correct by strand and frame info
    if strand == "-":
        seq = seq.reverse_complement()
        frames.reverse()
    seq = seq[frames[0]:]
    
    #get seq records
    #SeqRecord(Seq("MKQH", IUPAC.protein), id="YP_025292.1", name="HokC", description="toxic membrane protein")
    #function = "%snt %saa %s %s %s" % (len(seq), len(seq)/3.0, strand, str(boundaries), str(frames))
    cds = SeqRecord(seq, id=gene, name="", description=function)
    pep = SeqRecord(seq.translate(table=codonTable, to_stop=False), id=gene, name="", description=function)
    return cds, pep

def gff2fasta(gff, fasta, entireGene, codonTable, verbose):
    """Report gene, cds, peptide.
    not reporting 1000bp flanking intergenic sequence.
    """
    #load genome
    chr2seq = genome2dict(fasta)
    #load gff
    gene2position, contig2gene = load_gtf(gff)
    #get out streams
    genout = open(gff.name+".gene.fa", "w")
    cdsout = open(gff.name+".cds.fa", "w")
    pepout = open(gff.name+".pep.fa", "w")
    #process entries
    i = 0
    genes = set()
    for ci, contig in enumerate(sorted(contig2gene), 1):
        sys.stderr.write(" %s %s            \r" % (ci, contig))
        for s, e, feature, gene in contig2gene[contig]:
            i += 1
            contig, boundaries, strand, function, frames = gene2position[gene]
            #store CDS and peptide    
            cds, pep = get_cds_pep(chr2seq[contig], gene, boundaries, strand, function, frames, codonTable)
            cdsout.write(cds.format('fasta'))
            pepout.write(pep.format('fasta'))
            #get geneid and store gene
            geneid = ".".join(gene.split(".")[:-1])#; print geneid
            if geneid not in genes:
                genes.add(geneid)
                seq = chr2seq[contig][s-1:e]
                if strand == "-":
                    seq = seq.reverse_complement()
                gen = SeqRecord(seq, id=geneid, name="", description=function)
                genout.write(gen.format('fasta'))            
  
def main():
    
    usage   = "%(prog)s [options] -v" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v","--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-f", "--fasta",      type=argparse.FileType('r'),
                        help="genome stream")
    parser.add_argument("-g", "--gff",        type=argparse.FileType('r'),
                        help="annotation stream")
    parser.add_argument("--entireGene",       default=False, action="store_true", 
                        help="only full genes (start+stop codon)")
    parser.add_argument("-c", "--codonTable", default=12,
                        help="define codonTable  [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    gff2fasta(o.gff, o.fasta, o.entireGene, o.codonTable, o.verbose)
              
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
