#!/usr/bin/env python

import sys

out = sys.stdout
n = 4
pchrom = ''
i = chromi = pe = 0
for l in sys.stdin:
  chrom, s, e = l[:-1].split('\t')
  chrom = chrom.split()[0].split('|')[-1].split('.')[0]
  s, e = int(s), int(e)
  if not pchrom or chrom!=pchrom: 
    chromi = 1
    pchrom = chrom
    i = pe = 0
  out.write('%s.%s\t%s\t%s\n'%(chrom, chromi, s-pe, e-pe))
  i += 1
  if i>=n:
    i = 0
    chromi += 1
    pe = e
    
