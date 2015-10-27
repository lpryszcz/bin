#!/usr/bin/env python
# Increase intervals from BED by N nucleotides, strand specific. 

# USAGE: cat bed | bed2region.py 100000 > enlarged.bed

import os, sys

n = int(sys.argv[1])

for l in sys.stdin:
  ldata = l[:-1].split('\t')
  s, e = map(int, ldata[1:3])
  strand = ldata[5]
  #print ldata, (s,e,strand)
  if strand == "+":
    e += n
    ldata[2] = str(e)
  elif strand == "-":
    s -= n
    if s<0:
      s = 0
    ldata[1] = str(s)
  
  sys.stdout.write("\t".join(ldata)+"\n")

