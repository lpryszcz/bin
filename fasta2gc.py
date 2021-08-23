#!/usr/bin/env python2
# Return %GC of given sequence in FastA

# USAGE: cat fasta | fasta2gc.py 

import sys
from Bio import SeqIO
from Bio.SeqUtils import GC

for r in SeqIO.parse(sys.stdin, 'fasta'):
  print GC(r.seq)

