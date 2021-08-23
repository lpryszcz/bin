#!/usr/bin/env python2
###
# convert transcript positions into genomic positions. 
#
# USAGE: cat transcripts.txt | python transcript2genomic.py gtf > genomic.bed
# zcat GSE53693_allinfo_final_ORFs.txt.gz | cut -f3,4,9 | transcript2genomic.py ../../ref/DANRE.gtf > GSE53693.bed 

import sys
gtf = sys.argv[1]

def load_transcripts(gtf):
  """Return transcript dictonary"""
  tid2pos = {}
  for l in open(gtf):
    if l.startswith('#'):
      continue
    # 4       ensembl exon    52002   52120   .       -       .       gene_id "ENSDARG00000104632"; gene_version "1"; transcript_id "ENSDART00000166186"
    chrom, source, ftype, s, e, score, strand, frame, comment = l[:-1].split('\t')
    s, e = int(s)-1, int(e)
    if ftype != "exon":
      continue
    # get tid
    tid = filter(lambda x: x.startswith("transcript_id"), comment.split('; '))[0].split('"')[1]
    # store
    if tid not in tid2pos:
      tid2pos[tid] = []
    tid2pos[tid].append((chrom, s, e, strand))
  # sort
  for tid in tid2pos:
    tid2pos[tid].sort()
    if tid2pos[tid][0][-1] == "-":
      tid2pos[tid].reverse()
  return tid2pos
  
def get_genomic_coordinate(pos, tdata, length=1):
  """Return genomic coordinate for give transcript position"""
  # get transcript offset
  offset = 0
  for chrom, s, e, strand in tdata:
    elen = e - s
    if offset + elen > pos: 
      break
    offset += elen
  # get + / - coordinate
  if strand == "+":
    gs = s + pos - offset - 1
  else:
    gs = e - pos + offset
  return chrom, gs, gs+length, strand  

# load gtf
sys.stderr.write("Loading GTF...\n")
tid2pos = load_transcripts(gtf)
sys.stderr.write(" %s transcripts loaded!\n"%len(tid2pos))

sys.stderr.write("Parsing...\n")
missing = set()
for i, l in enumerate(sys.stdin, 1):
  if not i % 1000:
    sys.stderr.write(" %i   \r"%i)
  # unload line
  tid, pos, score = l[:-1].split('\t')[:3]
  try:
    pos = int(pos)
  except:
    sys.stderr.write("[WARNING] Cannot unload position: %s\n" % pos)
    continue
  # capture missing transcripts
  if tid not in tid2pos:
    if tid not in missing:
      #sys.stderr.write("[WARNING] Missing transcript: %s\n"%tid)
      missing.add(tid)
    continue
  # get genomic coordinate
  chrom, gs, ge, strand = get_genomic_coordinate(pos, tid2pos[tid])
  # report bed
  sys.stdout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom, gs, ge, tid, score, strand))
  
sys.stderr.write("%s transcripts were missing!\n"%len(missing))

