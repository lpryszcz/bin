#!/usr/bin/env python
"""GTF2BED
Report only transcripts
"""

import os, sys

def comment2transcript(comment):
  """Return gene id and gene name from comment"""
  k2v={}
  for kv in comment.strip(';').split(';'):
    kv.strip()
    k, v = kv.split(' "')
    k = k.strip()
    v = v.strip('"')
    k2v[k] = v
  #print k2v
  tid = k2v["transcript_id"]
  if "gene_name" in k2v:
    gene_name = k2v["gene_name"]
  else:
    gene_name = "-"
  return tid, gene_name


i=0
for line in sys.stdin:
  i+=1
  if line.startswith('#'):
    continue
  
  try:  
    contig,pred,f,s,e,dot,strand,score,comment = line.split('\t')
  except:
    sys.stderr.write( "Warning: Wrong line %s: %s\n" % ( i,str(line.split('\t')) ) )
    continue
    
  if f == 'transcript':
    #get geneID - gene_id "CPAG_00001"; transcript_id "CPAG_00001.1" -> CPAG_00001;
    #g = comment.split('"')[1]
    #t = comment.split('"')[3]
    tid, gene_name = comment2transcript(comment)
    s,e = int(s),int(e)
    #BED is 0-based, half-open
    s -= 1    
    print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig, s, e ,tid, score, strand, gene_name)
