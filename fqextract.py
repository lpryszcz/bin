#!/usr/bin/env python3
# Extract reads from fastq file that are present in list.
# USAGE: zcat reads.fq.gz | fqextract.py <(grep TRUE sequencing_summary.txt |cut -f2) | gzip > pass.fq.gz

import sys

out = sys.stdout
flag = False
ids = set(x.split()[0] for x in open(sys.argv[1], 'rt'))
sys.stderr.write(" %s read names loaded.\n"%len(ids))
k=0
for lineno, l in enumerate(sys.stdin, 1):
    if not lineno%1e5:
        sys.stderr.write(" %i processed with %i matche(s)\r"%(lineno/4, k/4))
    if lineno%4 == 1:
        flag = (l[1:].split()[0] in ids)
    if flag:
        out.write(l)
        k+=1
sys.stderr.write("%s extracted from %s reads\n"%(k/4, lineno/4))
