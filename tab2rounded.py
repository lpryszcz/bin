#!/usr/bin/env python
# Round and normalise expression values

# USAGE: cat txt | tab2rounded.py [1.0] > rounded.txt

import os, sys

norm = 1.0
if len(sys.argv)>1:
    norm = float(sys.argv[1])

ids = set()
for i, l in enumerate(sys.stdin):
    if not i:
        sys.stdout.write(l)
        continue
    ldata = l[:-1].split('\t')
    gid = ldata[0].replace(' ', '_')
    if gid in ids:
        sys.stderr.write('[WARNING] Duplicated ID: %s\n'%gid)
        #'''
        if gid.split('.')[-1].isdigit():
            gid = ".".join(gid.split('.')[:-1])+str(int(gid.split('.')[-1])+1)
        else:
            gid = gid+".1"
        #'''
    ids.add(gid)
    ndata = [gid,] + [str(int(round(float(x)*norm))) for x in ldata[1:]]
    sys.stdout.write("\t".join(ndata)+"\n")    