#!/usr/bin/env python
"""Selecte longest isoform. Uses augustus annotation (.t1/.t2)

USAGE: cat proteins.fa | isoforms2longest.py > proteins.longest.fa
"""

import sys
from Bio import SeqIO

#load fastas
id2seq = SeqIO.to_dict(SeqIO.parse(sys.stdin, 'fasta'))
sys.stderr.write("%s entries loaded\n"%len(id2seq))
#print id2seq
#get isoforms
gene2isoforms = {}
for tid in id2seq:
	gname = ".t".join(tid.split('.t')[:-1])
	if gname not in gene2isoforms:
		gene2isoforms[gname] = []
	gene2isoforms[gname].append(tid)
sys.stderr.write("%s genes recognised\nReporting longest isoforms...\n"%len(gene2isoforms))	
	
#return longest
for gname, tids in gene2isoforms.iteritems():
	longest = sorted(tids, key=lambda x: len(id2seq[x]), reverse=1)[0]
	sys.stdout.write(id2seq[longest].format('fasta'))

