#!/usr/bin/env python2
# Convert sequence to reverse complement. Assumes FastA input / output. 

# USAGE: cat fastq | fastx2reverse.py [fastq]

import sys
#from Bio import SeqIO

seqformat = "fasta"
if len(sys.argv)>1:
  seqformat = sys.argv[1]

def process_fasta(handle, out=sys.stdout):
  """Reverse seq in FastA format"""
  for l in sys.stdin:
    if l.startswith(">"):
      out.write(l)
    else:
      out.write("".join(reversed(l[:-1]))+"\n")
      
def process_fastq(handle, out=sys.stdout):
  """Reverse seq in FastQ format assuming 4-lines per entry"""
  for i, l in enumerate(sys.stdin):
    # seq and qual are odd lines
    if i%2:
      out.write("".join(reversed(l[:-1]))+"\n")
    # header and spacer (+) are even lines
    else:
      out.write(l)

if seqformat=="fasta":
  process_fasta(sys.stdin)
elif seqformat=="fastq":
  process_fastq(sys.stdin)
else:
  sys.exit("Sequence format (%s) not recognised!"%seqformat)
