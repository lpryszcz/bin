#!/usr/bin/env python
# convert exonerate gff into valid gtf

#USAGE: cat exonerate.gff | exonerategff2gtf.py > gtf

import sys
out = sys.stdout

def comment2dict(comment):
  """Unload comment line"""
  k2v = {}
  comment = comment.strip()
  if not comment:
    return k2v
  for kv in comment.split(';'):
    kv = kv.strip()
    if not kv:
      continue
    k = kv.split()[0]
    v = " ".join(kv.split()[1:])
    k2v[k] = v
  return k2v

gene_id = ""
genes = set()
for l in sys.stdin:
  ldata = l.split('\t')
  if l.startswith(('#',"Command line","Hostname","vulgar",">")) or len(ldata)!=9 or ldata[2] not in ('gene', 'cds'): #'exon'
    continue
  k2v = comment2dict(ldata[8])
  if "sequence" in k2v:
    gene_id = k2v["sequence"]
    i=1
    while "%s.%s"%(gene_id, i) in genes:
      i += 1
    gene_id = "%s.%s"%(gene_id, i)
    genes.add(gene_id)
  comment = 'gene_id "%s"; transcript_id "%s"'%(gene_id, gene_id)
  out.write("%s\t%s\n"%("\t".join(ldata[:8]), comment))
  
