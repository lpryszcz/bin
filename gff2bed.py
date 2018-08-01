#!/usr/bin/env python
"""GFF2BED"""

import os, sys

f2g = {}
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
    
  #add features
  if f not in f2g:
    f2g[f] = {}
    
  #if f=='gene':
  #ID=GENENAME;...
  g = comment.split('ID=')[1].split(";")[0]+"|"+f
  s,e = int(s),int(e)
  #BED is 0-based, half-open
  s -= 1
  print "%s\t%s\t%s\t%s\t%s\t%s" % ( contig,s,e,g,score,strand )
  
