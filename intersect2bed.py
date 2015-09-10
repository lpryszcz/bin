#!/usr/bin/env python
# Decrease intervals from BED by N nucleotides, strand specific. 
# and collapse overlapping intervals

# USAGE: bedtools intersect -a bed -b gtf | intersect2bed.py > collapsed.bed

import os, sys


def comment2gene(comment):
  """Return gene id and gene name from comment"""
  k2v={}
  for kv in comment.strip(';').split(';'):
    kv.strip()
    k, v = kv.split(' "')
    k = k.strip()
    v = v.strip('"')
    k2v[k] = v
  #print k2v
  gene_id = k2v["gene_id"]
  if "gene_name" in k2v:
    gene_name = k2v["gene_name"]
  else:
    gene_name = "-"
  return gene_id, gene_name

pldata = []
ids, names = [], []
for l in sys.stdin:
  ldata = l[:-1].split('\t')#; print ldata
  motif, score, strand, ref = ldata[3:7]
  # skip motifs without genic overlap
  if ref==".": 
    pldata = ldata
    continue
  #print ldata, (s,e,strand)
  # report collapsed matches
  if pldata and ldata[:4] != pldata[:4]:
    s, e = map(int, pldata[1:3])
    if strand == "+":
      e = s+len(motif)
      pldata[2] = str(e)
    elif strand == "-":
      s = e-len(motif)
      pldata[1] = str(s)
    sys.stdout.write("\t".join(pldata[:6]+["; ".join(ids), "; ".join(names)])+"\n")
    ids, names = [], []
  # get and store gene info    
  comment = ldata[14]
  gene_id, gene_name = comment2gene(comment)
  ids.append(gene_id)
  names.append(gene_name)
  
  pldata = ldata

if pldata and ldata[:4] != pldata[:4]:
  if strand == "+":
    e = s+len(motif)
    ldata[2] = str(e)
  elif strand == "-":
    s = e-len(motif)
    ldata[1] = str(s)
  sys.stdout.write("\t".join(pldata[:6]+["; ".join(ids), "; ".join(names)])+"\n")
  ids, names = [], []

