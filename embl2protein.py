#!/usr/bin/env python2
"""
Convert genbank to multifasta of proteins.

USAGE:
cat file.embl | embl2protein.py > file.faa

NOTE:
It's designed to work with embl files coming from EBI. gene is used as gene_id and transcript_id (locus_tag if gene not present).
Only entries having types in allowedTypes = ['gene','CDS','tRNA','tmRNA','rRNA','ncRNA'] are stored in GTF. Need to include exon processing.
No frame info is processed. Need to be included in order to process genes having introns!

AUTHOR:
Leszek Pryszcz
lpryszcz@crg.eu

Version 0.1

"""

import os, sys
from datetime import datetime
from Bio      import SeqIO

def embl2protein( source='embl2gtf',allowedTypes=set(['gene','CDS','tRNA','tmRNA','rRNA','ncRNA']) ):
  """
  """
  handle = sys.stdin
  for gb in SeqIO.parse( handle,'embl' ):
    acc     = gb.id #gb.name #gb.description # # 
    skipped = 0
    skippedTypes = set()
    for f in gb.features:
    
      #process only gene and CDS entries
      if f.type != 'CDS': #not in allowedTypes:
        #skipped += 1
        #skippedTypes.add( f.type )
        continue
      
      #generate comments field
      if 'locus_tag' in f.qualifiers:
        #use locul tag as gene_id/transcript_id
        gene_id = transcript_id = f.qualifiers['locus_tag'][0]
      else:
        sys.stderr.write( "Error: Neither `gene` nor `locus_tag` found for entry: %s\n" % '; '.join( str(f).split('\n') ) )
        continue
      
      comments = 'gene_id "%s"; transcript_id "%s"' % ( gene_id,transcript_id )
        
      if 'gene' in f.qualifiers:
        comments += '; gene_id "%s"' % f.qualifiers['gene'][0]
      if 'protein_id' in f.qualifiers:
        comments += '; protein_id "%s"' % f.qualifiers['protein_id'][0]
      
      #add external IDs  
      if 'db_xref' in f.qualifiers:
        for extData in f.qualifiers['db_xref']:
          comments += '; db_xref "%s"' % extData
      
      #code strand as +/- (in genbank 1 or -1)
      if int(f.strand)>0: strand = '+'
      else:               strand = '-'
      
      fasta = ">%s\n%s\n" % ( gene_id,f.extract(gb.seq).translate() )
      sys.stdout.write( fasta )
      
    #sys.stderr.write( "%s\tSkipped %s entries having types: %s.\n" % ( gb.id,skipped,', '.join(skippedTypes) ) )

if __name__=='__main__': 
  t0=datetime.now()
  embl2protein()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
