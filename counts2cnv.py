#!/usr/bin/env python
desc   = "Identify Copy Number Variants (CNVs)."
epilog = """Takes two files containing gene and read counts info (bam2counts.py output).
Note, gene counts has to be already normalised (RPKMs)! Use bam2counts.py for that;)

Author:
l.p.pryszcz@gmail.com

Barcelona, 25/04/2012
"""

import os, sys
from optparse import OptionParser
from datetime import datetime
from math     import log
from genome_annotation import load_gtf,load_gff
from snp2annotate import load_pfam, load_fasta_headers

def load_counts( fn,gene2counts={} ):
    """Return dictionary of gene counts."""
    firstIter = True
    if gene2counts:
        firstIter = False
    for l in open( fn ):
        l = l.strip()
        if not l or l.startswith('#'):
            continue
        lData = l.split("\t")
        if len(lData) > 2:
            contig,strand,start,stop,exon,gene,count,length,coverage,rpkm,log2rpkm1,gc = lData
            count      = rpkm
            coordinate = "%s:%s-%s%s" % (contig,start,stop,strand)
        else:
            gene,count = l.Data
            coordinate = ""

        try:
            count      = float( count )
        except:
            continue

        if firstIter:
            gene2counts[gene] = [coordinate,[]]
            
        gene2counts[gene][1].append( count )
    return gene2counts

def process( fnames,faa,pfam,gtf,log2th,splitFn,skipExons,verbose ):
    """main function
    """
    #load function annotation
    trans2ann = trans2pfam = {}
    if faa:
        trans2ann = load_fasta_headers( faa )
    if pfam:
        trans2pfam = load_pfam( pfam )

    ctg2cds,id2gene = {},{}
    if gtf:
        # load gtf/gff
        if gtf.endswith(".gff"):
            id2gene,ctg2cds = load_gff( gtf )
        else:
            id2gene,ctg2cds = load_gtf( gtf )
        if verbose:
            sys.stderr.write( "Loaded annotation of %s CDS from %s\n" % ( len(id2gene),o.gtf ) )
                
    #get samples names
    samples = []
    for fn in fnames:
        if splitFn:
            fn = fn.split(".")[0]        
        samples.append( fn )

    #load gene counts
    if verbose:
        sys.stderr.write( "Loading gene counts...\n" )
    gene2counts = {}
    for fn in fnames:
        if verbose:
            sys.stderr.write( " %s      \r" % fn )
        gene2counts = load_counts( fn,gene2counts )

    ## print results
    # header
    if verbose:
        sys.stderr.write( "Calculating...\n" )
    header = "#gene\tcoordinate\t%s" % samples[0]
    for s in samples[1:]:
        header += "\t%s\tlog2(%s/%s)" % ( s,s,samples[0] )
    header += "\tannotation\tpfam"
    print header
    
    # per gene scores
    for gene in sorted( gene2counts.keys() ):
        #if genes only requested then skip
        if skipExons:
            #check if exon, and skip if so
            if gene.split(".")[-1].isdigit():
                continue
        #
        coord,counts = gene2counts[gene]
        passed = False
        line   = "%s\t%s\t%.2f" % ( gene,coord,counts[0] )
        for c in counts[1:]:
            line += "\t%.2f" % c
            #ref 0
            if not counts[0]:
                line  += "\t+NA"
                passed = True
            elif not c:
                line  += "\t-NA"
                passed = True
            else:
                log2  = log( c*1.0/counts[0],2 )
                line += "\t%.2f" % log2
                #filter lines that contain log2 > than log2th or log2 < -log2th
                if log2th:
                    if not -log2th < log2 < log2th:
                        passed = True
                else:
                    passed = True
                
        #print only if passed filtering
        if passed:
            ann = pfam = ""
            if gene in id2gene:
                ann  = id2gene[gene][-1] #contig,cdsList,strand,function
            if gene in trans2ann:
                ann  = trans2ann[gene]
            if gene in trans2pfam:
                pfam = trans2pfam[gene]
            line += "\t%s\t%s" % ( ann,pfam )
            print line
    
def main():
    usage  = "%prog [options]"
    parser = OptionParser( usage=usage,description=desc,epilog=epilog,version="%prog 1.0" ) #allow_interspersed_args=True

    parser.add_option("-g", dest="gtf",
                      help="genome annotation gtf/gff [required if not -i]" )    
    parser.add_option("-l", dest="log2th",   default=0.0, type=float, 
                      help="report genes passing log2 cutoff [%default]")
    parser.add_option("-p", dest="pfam", default="", 
                      help="pfam tblout file") 
    parser.add_option("-q", dest="faa", default="", 
                      help="proteome fasta (to get protein annotation)")
    parser.add_option("-e", dest="skipExons",  default=False, action="store_true", 
                      help="report only gene loss [%default]")
    parser.add_option("-s", dest="splitFn",  default=False, action="store_true", 
                      help="split fname (sheet name) by dot")
    parser.add_option("-v", dest="verbose",  default=False, action="store_true" )
  
    ( o, fnames ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nFiles: %s\n" % ( o,fnames ) )

    for fn in fnames:
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )

    process( fnames,o.faa,o.pfam,o.gtf,o.log2th,o.splitFn,o.skipExons,o.verbose )
    
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
