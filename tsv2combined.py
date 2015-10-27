#!/usr/bin/env python
# Combine GO terms for each gene. 
# USAGE: cat biomart.txt | python tsv2combined.py > ensembl80.goterms.txt

import sys

gene2go = {}
for l in sys.stdin:
  if l.startswith("Ensembl"): #len(l.split('\t'))!=2:
    continue
  # read 
  gene, go = l[:-1].split('\t')
  # store gene
  if gene not in gene2go:
    gene2go[gene] = set()
  # store go
  if go: 
    gene2go[gene].add(go) 
  
# report
out = open("population.txt", "w")
for gene, gos in gene2go.iteritems():
  out.write("%s\n"%gene)
  if gos:
    sys.stdout.write("%s\t%s\n"%(gene, ";".join(gos)))

sys.stderr.write("Gene names written to: population.txt\n")

