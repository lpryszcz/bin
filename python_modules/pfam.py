#!/usr/bin/env python
"""This file store python modules to handle 
genome annotation formats like GTF, GFF, EMBL, GenBank.
hmmer tblout etc
"""
import sys

def load_pfam2go( fpath ):
    """Load pfam2go entries.
    Association file can be downloaded from: http://www.geneontology.org/external2go/pfam2go
    """
    pfam2go = {}
    for line in open(fpath):
        if line.startswith('!'):
            continue
        line  = line.strip()
        #Pfam:PF00001 7tm_1 > GO:G-protein coupled receptor protein signaling pathway ; GO:0007186
        pfam  = line.split()[0].split(':')[1]
        go    = line.split()[-1]
        #add entry
        if not pfam in pfam2go:
            pfam2go[pfam]=set()
        pfam2go[pfam].add( go )
      
    return pfam2go

def load_broad_pfam( gene2pfamFn,eTh=1e-05 ):
    """Return dictionary of list with GO mappings for each protein"""
    #PROTEIN_NAME  LOCUS  GENE_CONTIG  PFAM_ACC  PFAM_NAME  PFAM_DESCRIPTION  PFAM_START  PFAM_STOP  LENGTH  PFAM_SCORE  PFAM_EXPECTED
    prot2pfam = {}
    for line in open( gene2pfamFn ): 
        prot,gene,contig,pfam,pfam_name,desc,start,stop,length,score,e=line.split('\t')
        #check E-value
        e,score = float(e),float(score)
        if e>eTh:
            continue
        data=(pfam_name,desc,start,stop,length,score,e)
        if not prot in prot2pfam:
            prot2pfam[prot]={}
        prot2pfam[prot][pfam]=data
    return prot2pfam

def load_pfam_tblout( gene2pfamFn,eTh=1e-05 ):
    """Return dictionary of list with GO mappings for each protein."""
    #pfam_name,pfam,prot,gene,e,score,bias,
    prot2pfam = {}
    start=stop=length=0
    for line in open( gene2pfamFn ):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        #it's multiple spaces table
        data=line.split()
        try:
            pfam_name,pfam,prot,gene,e,score,bias = data[:7]
        except:
            sys.stderr.write( "Error: load_pfam_tblout: Cannot parse line: %s\n" % data )
        desc = " ".join( data[18:] )
        #check E-value
        e,score = float(e),float(score)
        if e>eTh:
            continue
        #correct phylomedb isoforms Phy003WV3G_9999994_2 > Phy003WV3G_9999994
        if prot.startswith("Phy") and prot.count("_")>1:
            prot = "_".join(prot.split("_")[:2])
        data=(pfam_name,desc,start,stop,length,score,e)#; print data
        if not prot in prot2pfam:
            prot2pfam[prot]={}
        prot2pfam[prot][pfam]=data
    return prot2pfam
    

    