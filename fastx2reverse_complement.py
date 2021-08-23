#!/usr/bin/env python2
# Convert sequence to reverse complement. Assumes FastA input / output. 

# USAGE: cat fastq | fastx2reverse_complement.py [fastq]

import sys
from Bio import SeqIO

out = sys.stdout

seqformat = "fasta"
if len(sys.argv)>1:
  seqformat = sys.argv[1]

for r in SeqIO.parse(sys.stdin, seqformat):
  rc = r.reverse_complement()
  rc.id, rc.description = r.id, r.description
  out.write(rc.format(seqformat))

