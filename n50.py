#!/usr/bin/env python
# Calculate N50 of given int

# USAGE: cat int.txt | n50.py [50]

import sys
#import numpy as np

def get_n50(handle, n=50):
  """Return N50 statistics"""
  values = sorted([int(i) for i in handle], reverse=1)
  nsize  = sum(values) * n / 100.0
  cumulative = 0
  for v in values:
    cumulative += v
    if cumulative >= nsize:
      break
  print "N%s = %s"%(n, v)

if __name__ == "__main__":
  if len(sys.argv)>1:
    n = int(sys.argv[1])
  get_n50(sys.stdin, n)
