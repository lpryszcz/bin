#!/usr/bin/env python
# Report bed with contigs to chromosomes associations. Store only best!
# USAGE: 
#  lastal -P 4 chromosomes.fa contigs.fa $f | last-split - | maf-convert tab - | ./tab2chromosome.py > contigs.bed
#  awk '{chr+=$3-$2; score+=$5}END{print 100*score/chr"%", chr, score}' contigs.bed # to get the identity

import sys

out = sys.stdout
q2chr = {}
minBestScoreFrac = 0.9

for l in sys.stdin:
  if l.startswith('#'): continue
  score, t, tstart, tlen, tstrand, tsize, q, qstart, qlen, qstrand, qsize = l.split('\t')[:11]
  score, tstart, tlen, qstart, qlen, qsize = map(int, (score, tstart, tlen, qstart, qlen, qsize))
  # update start if reverse - this will provide q coords in forward strand
  if qstrand == "-":
    qstart = qsize-qstart-qlen
  if q in q2chr:
    if score > q2chr[q][-1]: #< minBestScoreFrac*q2chr[q][0][-1]: 
      q2chr[q] = (t, tstart, tstart+tlen, qstrand, score)
  else:
    q2chr[q] = (t, tstart, tstart+tlen, qstrand, score)
    #q2chrs[q].append((score, t, tstart, tstart+tlen, qstrand))

for q, (t, start, end, qstrand, score) in sorted(q2chr.iteritems(), key=lambda x: x[1]):
  out.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(t, start, end, q, score, qstrand))

