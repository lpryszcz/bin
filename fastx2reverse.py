#!/usr/bin/env python
# Convert sequence to reverse complement. Assumes FastA input / output. 

# USAGE: cat fastq | fastx2reverse.py [fastq]

import sys
#from Bio import SeqIO

out = sys.stdout

seqformat = "fasta"
if len(sys.argv)>1:
  seqformat = sys.argv[1]
if seqformat != "fasta":
  sys.exit("No FastQ support")

for l in sys.stdin:
  if l.startswith(">"):
    out.write(l)
  else:
    out.write("".join(reversed(l[:-1]))+"\n")


