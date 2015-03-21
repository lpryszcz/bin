#!/usr/bin/env python
# Extract reads from fastq file.

import sys 
ids, lineno, flag = set('@' + x for x in open(sys.argv[1])), 0, False
for l in sys.stdin:
    lineno += 1
    if lineno%4 == 1: flag = (l not in ids)
    if flag: print l,
