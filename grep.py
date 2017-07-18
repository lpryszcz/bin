#!/usr/bin/env python
### Fast `grep -f` implementation for GTF files
# USAGE: grep.py file_with_elements file_to_grep

import sys

if len(sys.argv)<3:
  sys.exit("Provide two arguments!")

fn1, fn2 = sys.argv[1:3]

if fn2 == "-":
  handle2 = sys.stdin
else:
  handle2 = open(fn2)
  
lines = set(l[:-1] for l in open(fn1))

for i, l in enumerate(handle2, 1):
  if not i%100:
    sys.stderr.write(" %i      \r"%i)
  if lines.intersection([e.strip(';').strip('"') for e in l.split(' ',)]):
    sys.stdout.write(l)
    

