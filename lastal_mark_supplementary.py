#!/usr/bin/env python3

import sys

handle = sys.stdin
out = sys.stdout

def report_read(out, algs):
    """Report read marking supplementary alignments"""
    # get sorted algs for the same read, starting with the highest mapq
    # if more than 1 has the same mapq, get the one with longer alg length as primary
    algs = sorted(algs, key=lambda a: (int(a[4]), len(a[9])), reverse=True)
    for i, a in enumerate(algs):
        # mark supplementary algs - simply not the primary
        if i:
            a[1] = str(int(a[1])+2048)
        # and report
        out.write("\t".join(a))

algs = []
for l in handle:
    if l.startswith("@"):
        out.write(l)
        continue
    ldata = l.split('\t')
    if algs and ldata[0]!=algs[-1][0]:
        report_read(out, algs)
        algs = []
    algs.append(ldata)

if algs:
    report_read(out, algs)
    
