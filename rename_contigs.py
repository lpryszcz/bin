#!/usr/bin/env python
# Unify the fasta headers. Solve problems with | and other non-standard characters that cause some programs to fail. 
# USAGE: cat contigs.fa | rename_contigs.py [scaffold] > renamed.fa


import os,sys
from datetime import datetime
from Bio import SeqIO

def main(basename="cotnig"):
  
  for i, r in enumerate(SeqIO.parse(sys.stdin, 'fasta'), 1):
  	contig = "%s%5i" % (basename, i)
  	contig = contig.replace(' ', '0')
  	sys.stdout.write(">%s\n%s\n" % (contig, "\n".join(str(r.seq[j:j+60]) for j in range(0, len(r), 60))))

if __name__=='__main__': 
  t0=datetime.now()
  basename = "contig"
  if len(sys.argv)>1:
    basename = sys.argv[1]
  main(basename)
  dt=datetime.now()-t0
  sys.stderr.write( "Time elapsed: %s\n" % dt )
