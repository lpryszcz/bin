#!/usr/bin/env python
### Fast `grep -f` implementation for GTF files. If -s, then IDs are split by `.`.
# USAGE: grep.py [-s] file_with_elements file_to_grep

import sys

params = filter(lambda x: x.startswith("-"), sys.argv[1:])
fnames = filter(lambda x: not x.startswith("-"), sys.argv[1:])

split = False
if "-s" in params:
  split = True
  
if len(fnames)!=2:
  sys.exit("Provide two files!")

fn1, fn2 = fnames[:2]

if fn2 == "-":
  handle2 = sys.stdin
else:
  handle2 = open(fn2)

def dotsplit(s, split=True):
  """Split string by dot and return first element"""
  if split:
    return s.split('.')[0]
  return s

lines = set(dotsplit(l[:-1] ,split) for l in open(fn1))

for i, l in enumerate(handle2, 1):
  if not i%100:
    sys.stderr.write(" %i      \r"%i)
  if lines.intersection([dotsplit(e.strip(';').strip('"'), split) for e in l.split(' ',)]):
    sys.stdout.write(l)
    

