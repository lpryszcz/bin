#!/usr/bin/env python2
# Extract reads from fastq file that are present in list.

import gzip, sys

rnameFn = sys.argv[1]
if rnameFn.endswith(".gz"):
    rnameFile = gzip.open(rnameFn)
else:
    rnameFile = open(rnameFn)

out = sys.stdout
ids, flag = set('@' + x.split()[0] for x in rnameFile), False
sys.stderr.write(" %s read names loaded.\n"%len(ids))
k=0
for lineno, l in enumerate(sys.stdin, 1):
    if not lineno%1e5:
        sys.stderr.write(" %i processed with %i matche(s)\r"%(lineno/4, k/4))
    if lineno%4 == 1:
        flag = (l.split()[0] in ids)
    if flag:
        out.write(l)
        k+=1
    
sys.stderr.write("%s extracted from %s reads\n"%(k/4, lineno/4))
