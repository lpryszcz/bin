#!/usr/bin/env python2
# Report .bed of softmasked sequences
# USAGE: cat genome.fa | fasta2softmasked.py > repeats.bed

import re, sys
from Bio import SeqIO

pat=re.compile('[actgn]+')
out = sys.stdout

for i, r in enumerate(SeqIO.parse(sys.stdin, 'fasta'), 1):
    sys.stderr.write(' %s %s %s bp        \r'%(i, r.id, len(r)))
    for m in pat.finditer(str(r.seq)):
        s, e = m.span()
        out.write("%s\t%s\t%s\n"%(r.id, s, e))
        
