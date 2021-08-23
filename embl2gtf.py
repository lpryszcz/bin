#!/usr/bin/env python2
"""
Convert embl to gtf.

USAGE:
cat in.embl | embl2gtf.py > out.gtf

NOTE:
It's designed to work with gb files coming from EMBL. gene is used as gene_id and transcript_id (locus_tag if gene not present).
Only entries having types in allowedTypes = ['gene','CDS','tRNA','tmRNA','rRNA','ncRNA'] are stored in GTF. Need to include exon processing.
No frame info is processed. Need to be included in order to process genes having introns!
In addition, some embl files misses gene entries - then gene is considered as entire CDS/*RNA

AUTHOR:
Leszek Pryszcz
lpryszcz@crg.eu

Version 0.1

"""

import os, sys, urllib
from datetime import datetime
from Bio      import SeqIO

def _get_unique_id( gene_id,products,i=1 ):
  """ """
  new_id = '%s.%s' % ( gene_id,i )
  while new_id in products:
    i += 1
    new_id = '%s.%s' % ( gene_id,i )
  
  return new_id

def embl2gtf( source='embl2gtf',allowedTypes=set(['CDS', 'tRNA', 'tmRNA', 'rRNA', 'ncRNA']) ): #'mRNA', 'gene', 
  """
  """
  #make sure no duplicated products
  products = set()
  
  handle = sys.stdin
  for r in SeqIO.parse( handle,'embl' ):
    acc     = r.name # r.id r.description
    skipped = 0
    skippedTypes = set()
    for f in r.features:
    
      #process only gene and CDS entries
      if f.type not in allowedTypes:
        skipped += 1
        skippedTypes.add( f.type )
        continue
      
      #generate comments field
      if 'locus_tag' in f.qualifiers:
        #use locul tag as gene_id/transcript_id
        gene_id = transcript_id = f.qualifiers['locus_tag'][0]
      elif 'gene' in f.qualifiers:
        #use gene name as gene_id/transcript_id
        gene_id = transcript_id = f.qualifiers['gene'][0]
      elif 'product' in f.qualifiers:
        #use locul tag as gene_id/transcript_id
        gene_id = transcript_id = f.qualifiers['product'][0]
      else:
        sys.stderr.write( "Error: Neither `gene` nor `locus_tag` found for entry: %s\n" % '; '.join( str(f).split('\n') ) )
        continue
      
      #make sure no duplicated entries in GTF  
      if gene_id in products:
        #sys.stderr.write( "Warning: Duplicated entry found: %s\n" % '; '.join( str(f).split('\n') ) )
        continue
        #gene_id = transcript_id = _get_unique_id( gene_id,products ) #continue
      
      products.add(gene_id)
      
      comments = 'gene_id "%s"; transcript_id "%s"' % ( gene_id,transcript_id )
      #if not 'systematic_id' in f.qualifiers:
      #  continue
      #comments = 'gene_id "%s"; transcript_id "%s"' % ( f.qualifiers['systematic_id'][0],f.qualifiers['systematic_id'][0] ) 
        
      if 'protein_id' in f.qualifiers:
        comments += '; protein_id "%s"' % f.qualifiers['protein_id'][0]
        
      if 'note' in f.qualifiers:
        comments += '; note "%s"' % f.qualifiers['note'][0]
        
      if 'product' in f.qualifiers:
        if f.qualifiers['product'][0] != 'hypothetical protein':
          comments += '; product "%s"' % urllib.quote(f.qualifiers['product'][0])
      
      #add external IDs  
      if 'db_xref' in f.qualifiers:
        for extData in f.qualifiers['db_xref']:
          comments += '; db_xref "%s"' % extData
      
      #code strand as +/- (in genbank 1 or -1)
      if int(f.strand)>0: strand = '+'
      else:               strand = '-'
      
      #define gb
      """
      seqname - The name of the sequence. Must be a chromosome or scaffold.
      source - The program that generated this feature.
      feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
      start - The starting position of the feature in the sequence. The first base is numbered 1.
      end - The ending position of the feature (inclusive).
      score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
      strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
      frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
      comments - gene_id "Em:U62317.C22.6.mRNA"; transcript_id "Em:U62317.C22.6.mRNA"; exon_number 1
      """
      s,e = f.location.start.position+1,f.location.end.position
      gtf = '%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s\n' % ( acc,source,'gene',s,e,strand,comments )
      #add start and stop codon
      if strand=='+':
        gtf += '%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s\n' % ( acc,source,'start_codon',s,s+3,strand,comments )
        gtf += '%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s\n' % ( acc,source,'stop_codon', e-3,e,strand,comments )
      else:
        gtf += '%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s\n' % ( acc,source,'stop_codon',s,s+3,strand,comments )
        gtf += '%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s\n' % ( acc,source,'start_codon', e-3,e,strand,comments )
      
      if f.sub_features:
        for subf in f.sub_features:
          gtf += '%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t%s\n' % ( acc,source,f.type,subf.location.start.position+1,subf.location.end.position,strand,comments )
      else:       #use entire entry as CDS
        gtf += '%s\t%s\t%s\t%s\t%s\t.\t%s\t0\t%s\n' % ( acc,source,f.type,f.location.start.position+1,f.location.end.position,strand,comments ) #f.frame,
      sys.stdout.write( gtf )
      
    sys.stderr.write( "%s\tSkipped %s entries having types: %s.\n" % ( r.name,skipped,', '.join(skippedTypes) ) )

if __name__=='__main__': 
  t0=datetime.now()
  embl2gtf()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
