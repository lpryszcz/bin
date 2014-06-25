#!/usr/bin/env python
#Report number of LOH over 100bp and the longest ones.

import sys
from Bio import SeqIO

fn1, fn2 = sys.argv[1:3]

lohs = []
for r1, r2 in zip(SeqIO.parse(fn1,'fasta'),SeqIO.parse(fn2,'fasta')):
  sys.stderr.write(' %s   \r'%r1.id)
  k = 0
  for b1, b2 in zip(r1, r2):
    if b1!=b2:
      if k: lohs.append(k)
      k = 0
    else:
      k += 1
  if k: lohs.append(k)
  
print "The longest LOH: %s bp" % max(lohs)
lohs100 = filter(lambda x: x>=100, lohs)
print "%s LOHs > 100bp: %s bp" %(len(lohs100),sum(lohs100))
print "%s LOHs: %s bp" %(len(lohs),sum(lohs))
