#!/usr/bin/env python
"""Report synteny stats from last tab file

Note, it may fail at inversions that are in separate blocks

USAGE: zcat CANORvsCDC317.tsv.gz | python tab2synteny.py

"""

import sys

minLen = 1000
overlapTh = 0.1 

def get_overlap(data1, data2, overlapTh):
  """Return False if two algs are not overlapping. 
  Return 1 if data1 is larger, 2 otherwise."""
  t1, ts1, tlen1 = data1[:3]
  t2, ts2, tlen2 = data2[:3]
  if t1==t2 and ts1+tlen1 - ts2 > overlapTh*max((tlen1, tlen2)):
    if tlen1>tlen2:
      return 1
    else:
      return 2

def get_last_size(sblocks):
  #end - start
  return sblocks[-1][-1][1]+sblocks[-1][-1][2] - sblocks[-1][0][1]

print "Parsing matches..."
i = k = 0
matches = []
seq2size = {}
for l in sys.stdin:
  if l.startswith('#'):
    continue
  i += 1
  name, t, ts, tlen, tstrand, tsize, q, qs, qlen, qstrand, qsize, algstring = l[:-1].split('\t')
  ts, tlen, tsize, qs, qlen, qsize = map(int, (ts, tlen, tsize, qs, qlen, qsize))
  #store seq size
  seq2size[ts] = tsize
  seq2size[qs] = qsize
  #skip if too short
  if tlen<minLen or qlen<minLen:
    continue
  #store 
  data = (t, ts, tlen, tstrand, q, qs, qlen, qstrand)
  #
  if matches and get_overlap(matches[-1], data, overlapTh):
    if get_overlap(matches[-1], data, overlapTh)==2:
      matches[-1] = data
  else:
    matches.append(data)

rAlg, qAlg = sum(m[2] for m in matches), sum(m[6] for m in matches)
print " %s / %s matches kept:\n  %s bp of ref aligned\n  %s bp of query aligned"%(len(matches), i, rAlg, qAlg)

print "Defining synteny blocks..."
minSlen, minSsize = 3, 10000
pdata = matches[0]
sblocks = [] #[[matches[0]]]
for data in matches:
  t, ts, tlen, tstrand, q, qs, qlen, qstrand = data
  #if different chromosome than the last synteny block
  if sblocks and t != sblocks[-1][-1][0] or sblocks and q!=sblocks[-1][-1][4]:
    #if synteny block is too short, remove it
    if len(sblocks[-1])<minSlen or get_last_size(sblocks) < minSsize:
      sblocks.pop(-1)
    #store current match as new synteny block 
    ###CHECK THE DISTANCE
    if sblocks and t != sblocks[-1][-1][0] or sblocks and q!=sblocks[-1][-1][4]:
      sblocks.append([data])
    else:
      if not sblocks:
        sblocks.append([])
      sblocks[-1].append(data)
  #extend last block
  else:
    if not sblocks:
      sblocks.append([])
    sblocks[-1].append(data)

#remove last block if to short
if sblocks and len(sblocks[-1])<minSlen \
   or sblocks and get_last_size(sblocks) < minSsize:
  sblocks.pop(-1)

print "#ID\tTarget\tQuery\tTarget alg length\tQuery alg length\tAligned fragments\tInversions"
Inversions = 0
tLength = 0
qLength = 0
for i, block in enumerate(sblocks, 1):
  tcoord="%s:%s-%s"%(block[0][0], block[0][1], block[-1][1]+block[-1][2])
  qcoord="%s:%s-%s"%(block[0][4], block[0][5], block[-1][5]+block[-1][6])
  tlength = block[-1][1]+block[-1][2] - block[0][1]
  tLength += tlength
  qlength = abs(block[-1][5]+block[-1][6] - block[0][5])
  qLength += qlength
  #count inversions
  strands = [block[0][-1]]
  for t, ts, tlen, tstrand, q, qs, qlen, qstrand in block[1:]:
    if qstrand != strands[-1]:
      strands.append(qstrand)
  inversions = len(strands)
  Inversions += inversions
  print "block_%s\t%s\t%s\t%s\t%s\t%s\t%s"%(i, tcoord, qcoord, tlength, qlength, len(block), inversions)
print "%s inversions found.\n  %s bp of ref in synteny blocks\n  %s bp of query in synteny blocks"%(Inversions, tLength, qLength)
