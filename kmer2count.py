#!/usr/bin/env python3
# Report counts of mers from FastA file
#USAGE: mer2count.py FastA [5 10]

import os, sys, pysam
#from Bio import SeqIO

n, m = 5, 10
fastafn = sys.argv[1]
if len(sys.argv)>3:
    n, m = map(int, sys.argv[-2:])

faidx = pysam.FastaFile(fastafn)
for k in range(n, m+1):
    mers = set()
    for ref in faidx.references:
        seq = faidx[ref]
        mers.update(set(seq[s:s+k] for s in range(len(seq)-k)))
    print(k, len(mers), 4**k, len(mers)/4**k)
