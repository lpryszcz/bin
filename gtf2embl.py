#!/usr/bin/env python
desc="""Report embl for gtf file.
Note: require Biopython 1.62+
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 14/06/2013
"""

import argparse, gzip, os, resource, sys
from datetime import datetime
from Bio          import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation, AfterPosition, BeforePosition, CompoundLocation
from genome_annotation import load_gtf

def get_locations(CDSs, start, end, strand):
    """Return mRNA and CDS locations
    CDS has exact boundaries, while mRNA not.
    """
    #gff is 1-based, gb also, but sf is 0-based
    if len(CDSs)>1:
        parts, mrnaparts = [], []
        for cdsi, (s, e) in enumerate(CDSs):
            parts.append(FeatureLocation(s-1, e, strand=strand))
            if cdsi==0:
                mrnaparts.append(FeatureLocation(BeforePosition(s-1), e, strand=strand))
            elif cdsi == len(CDSs)-1:
                mrnaparts.append(FeatureLocation(s-1, AfterPosition(e), strand=strand))
            else:
                mrnaparts.append(FeatureLocation(s-1, e, strand=strand))
        cdsloc  = CompoundLocation(parts)
        mrnaloc = CompoundLocation(parts)
    else:
        cdsloc  = FeatureLocation(start-1, end, strand=strand)
        mrnaloc = FeatureLocation(BeforePosition(start-1), AfterPosition(end), strand=strand)
    return cdsloc, mrnaloc

def gene2features(r, gene, gene2position, gene2product, start, end, gcode, partialyes, verbose):
    """
    """
    contig, CDSs, gffstrand, function, frames = gene2position[gene]
    if gffstrand in ('1','+'):
        strand = +1
    else:
        strand = -1
        CDSs.reverse()
    '''#add stop codon if not partial seq
    if strand==1 and CDSs[-1][1]+3 <= len(r.seq):
            CDSs[-1][1] += 3
    elif strand==-1 and CDSs[0][0]-3 > 0:
        CDSs[0][0] -= 3'''
    cdsloc, mrnaloc = get_locations(CDSs, start, end, strand)
    #add gene
    geneid = gene #".".join(gene.split('.')[:-1])
    #get product
    product = "hypothetical protein"
    if geneid in gene2product:
        product = gene2product[geneid]    
    if gene.endswith('.t1'):
        sf = SeqFeature(FeatureLocation(BeforePosition(start-1),AfterPosition(end)), strand=strand, type='gene', id=geneid)
        sf.qualifiers={"locus_tag": geneid, "gene": geneid, "product": product}
        r.features.append(sf)
    #get mRNA sf
    sf = SeqFeature(mrnaloc, type='mRNA', id=gene)
    sf.qualifiers={"locus_tag": geneid, "gene": geneid, "product": product} #"protein_id": gene
    r.features.append(sf)
    #get CDS sf
    sf = SeqFeature(cdsloc, type='CDS', id=gene)
    #get translation
    seq = sf.extract(r.seq)
    aa = str(seq.translate(table=gcode))
    #solve non-triplets issue
    if len(seq) % 3:
        if strand==1:
            end   -= len(seq) % 3
        else:
            start += len(seq) % 3
    ##check for partial sequence - no M as first or no * as last aa
    partial = 0
    #both ends partial
    if aa[0]!="M" and aa[-1]!="*":
        partial = 1
        sf.location = FeatureLocation(BeforePosition(start-1),AfterPosition(end))
    #left end partial
    elif aa[0]!="M" and strand==1 or aa[-1]!="*" and strand==-1:
        partial = 1                
        sf.location = FeatureLocation(BeforePosition(start-1),end)
    #right end partial
    elif aa[-1]!="*" and strand==1 or aa[0]!="M" and strand==-1:
        partial = 1
        sf.location = FeatureLocation(start-1,AfterPosition(end))
    #strip stop codon
    aa = aa.strip("*")
    #replace internal stop codons by X
    if "*" in aa:
        if verbose:
            sys.stderr.write("[Warning] Stop codon(s) in: %s. Skipped!\n" % gene)
        return r
        #aa = aa.replace("*","X")
    sf.qualifiers = {'transl_table': gcode, "locus_tag": geneid, "gene": geneid, "product": product, "translation": aa} #"protein_id": gene,
    if function:
        sf.qualifiers['note'] = function
    #inform about partial entries
    if partial:
        #skip if not partial are allowed
        if not partialyes:
            return r
        if aa[0]!="M":
            sf.qualifiers['codon_start'] = 1
        sf.qualifiers['product']    += ", partial cds"
        if verbose:
            sys.stderr.write("[Warning] Partial sequence: %s\n" % (gene,))
            #sys.stderr.write("[Warning] Partial sequence: %s %s\n" % (gene,sf))
    #add to features
    r.features.append(sf)
    return r 

def load_gene2product(handle, verbose=0):
    """Return gene2procuct dictionary"""
    if verbose:
        sys.stderr.write("Loading gene2product...\n")
    gene2product = {}
    for l in handle:
        l = l[:-1]
        gene, product = l.split('\t')
        gene2product[gene] = product
    if verbose:
        sys.stderr.write(" Products for %s genes loaded!\n" % len(gene2product))
    return gene2product
    
def gtf2embl(fasta, gtf, output, project, organism, strain, taxid, moltype, \
             topology, tech, gcode, productsfile, partial, verbose):
    """
    """
    if verbose:
        sys.stderr.write("Loading fasta...\n")
    #load fastas
    c2s = SeqIO.to_dict(SeqIO.parse(fasta,'fasta'))
    contigs = sorted(c2s, key=lambda x: len(c2s[x]), reverse=True)
    if verbose:
        sys.stderr.write(" %s fastas loaded: %s ...\n" % (len(contigs),", ".join(contigs[:5])))
    
    if verbose:
        sys.stderr.write("Loading annotation...\n")
    #load annotation
    gene2position, contig2gene = load_gtf(gtf.name, partial)
    if verbose:
        sys.stderr.write(" %s genes for %s fasta loaded.\n" % (len(gene2position), len(contig2gene)))

    gene2product = {}
    if productsfile:
        gene2product = load_gene2product(productsfile, verbose)
        
    #process contigs starting from the largets
    if verbose:
        sys.stderr.write("Generating embl file...\n")
    for i, c in enumerate(contigs, 1):
        if verbose:
            sys.stderr.write(" %s / %s  %s           \n" % (i, len(contigs), c))
        #set alphabet
        r = c2s[c]
        r.id = 'XXX'
        r.dbxrefs.append('Project:%s' % project)
        r.seq.alphabet = IUPAC.ambiguous_dna
        r.annotations  =  {'accessions': ['XXX'], "db_xref": "taxon:%s" %taxid, 
                           'organism': organism, 'sequence_version': 'XXX' }
        #'strain': strain, 'mol_type': moltype, "topology": topology }
        #add description
        """
    In [46]: embl.annotations
    Out[46]: 
    {'accessions': ['HE605202'],
     'data_file_division': 'FUN',
     'organism': 'Candida parapsilosis',
     'references': [Reference(title=';', ...),
      Reference(title='Transcriptional landscape of the pathogenic yeast Candida parapsilosis', ...)],
     'sequence_version': 1,
     'taxonomy': ['Eukaryota',
      'Fungi',
      'Dikarya',
      'Ascomycota',
      'Saccharomycotina',
      'Saccharomycetes',
      'Saccharomycetales',
      'mitosporic Saccharomycetales',
      'Candida']}

    SeqFeature(self, location=None, type='', location_operator='', strand=None,
            id='<unknown id>', qualifiers=None, sub_features=None, ref=None,
            ref_db=None)        
            """
        #add source features
        qualifiers = {'organism': organism, 'strain': strain, 'mol_type': moltype, 
                      'db_xref': 'taxon:%s' %taxid }
        sf = SeqFeature(location=FeatureLocation(0,len(r.seq)), strand=1, type='source', qualifiers=qualifiers) 
        r.features.append(sf)
        #check if any annotation
        if c not in contig2gene:
            if verbose:
                sys.stderr.write("[Warning] No annotation for: %s\n" % c)
            return
        elif verbose:
            sys.stderr.write("  %s genes predicted.\n" % len(contig2gene[c]))
        #add genes features
        for start, end, feature, gene in contig2gene[c]:
            r = gene2features(r, gene, gene2position, gene2product, start, end, gcode, partial, verbose)
        output.write(r.format('embl'))

def main():
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",   default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-f", dest="fasta",     required=True, type=file, 
                        help="fasta stream")
    parser.add_argument("-g", dest="gtf",       required=True, type=file, 
                        help="gtf stream")
    parser.add_argument("-o", dest="output",    default=sys.stdout, type=argparse.FileType("w"),
                        help="output stream     [stdout]")
    '''parser.add_argument("--functions", type=file, default=None, 
                        help="funtional annotations")    
    parser.add_argument("--goterms",   type=file, default=None, 
                        help="GO annotations")'''    
    parser.add_argument("--products",  type=file, default=None, 
                        help="products for genes") 
    parser.add_argument("--partial",   default=False, action="store_true",
                        help="allow partial genes") 
    sample = parser.add_argument_group('Sample information')
    sample.add_argument("-p", dest="project",   default='PRJEB1685', type=str,
                        help="organims          [%(default)s]")
    sample.add_argument("-r", dest="organism",  default='Candida parapsilosis', type=str,
                        help="organims          [%(default)s]")
    sample.add_argument("-s", dest="strain",    required=True, 
                        help="strain            [%(default)s]")
    sample.add_argument("-t", dest="taxid",     required=True, type=int, 
                        help="taxid")   
    sample.add_argument("-m", dest="moltype",   default='genomic DNA', 
                        help="molecule type     [%(default)s]")
    sample.add_argument("-w", dest="topology",  default='linear',
                        help="define topology   [%(default)s]")   
    sample.add_argument("-u", dest="tech",      default='wgs',
                        help="seq. technique    [%(default)s]")   
    sample.add_argument("-c", dest="gcode",     default=1, type=int,
                        help="genetic code      [%(default)s]") #transl_table=12
                         
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "[gtf2embl] Options: %s\n" % str(o) )

    gtf2embl(o.fasta, o.gtf, o.output, o.project, o.organism, o.strain, o.taxid, \
             o.moltype, o.topology, o.tech, o.gcode, o.products, o.partial, o.verbose)
	
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\n[gtf2embl] Ctrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write( "[gtf2embl] Time elapsed: %s\n" % dt )
