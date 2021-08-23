#!/usr/bin/env python2
"""GTF2BED"""

import os, sys

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
    
  if f == 'gene':
    #get geneID - gene_id "CPAG_00001"; transcript_id "CPAG_00001.1" -> CPAG_00001;
    g = comment.split('"')[1]
    t = comment.split('"')[3]
    s,e = int(s),int(e)
    #BED is 0-based, half-open
    s -= 1    
    print "%s\t%s\t%s\t%s\t%s\t%s" % ( contig,s,e,t,score,strand )
