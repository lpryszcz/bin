#!/usr/bin/env python
# Extract reads from fastq file that are present in list.

import gzip, sys

rnameFn = sys.argv[1]
if rnameFn.endswith(".gz"):
    rnameFile = gzip.open(rnameFn)
else:
    rnameFile = open(rnameFn)

ids, flag = set('@' + x.split()[0] for x in rnameFile), False
sys.stderr.write(" %s read names loaded.\n"%len(ids))
k=0
for lineno, l in enumerate(sys.stdin, 1):
    if lineno%1e5:
        sys.stderr.write(" %i\r"%lineno)
    if lineno%4 == 1:
        flag = (l.split()[0] in ids)
    if flag:
        print l,
        k+=1
    
sys.stderr.write("%s extracted from %s reads\n"%(k/4, lineno/4))
