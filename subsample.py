#!/usr/bin/env python2
# Report up to N random lines from input that are divided by last column
# USAGE: bedtools intersect -bed -F 0.6 -wb -a $f -b <(bedtools makewindows -g $ref.fai -w 1000 | grep chikungunya) | cut -f4,15 | subsample.py 1000 | cut -f1 | sort | uniq > $f.ids

import os, random, sys

N = int(sys.argv[1]) if len(sys.argv)>1 else 100#; print(N)
out = sys.stdout

cat2line = {}
for l in sys.stdin:
    cat = l[:-1].split("\t")[-1]
    if cat not in cat2line:
        cat2line[cat] = []
    cat2line[cat].append(l)

for cat, lines in cat2line.items():
    out.write("".join(random.choices(lines, k=N)))

