#!/usr/bin/env python
# Round expression values

# USAGE: cat txt | ./tab2rounded.py > rounded.txt

import os, sys

ids = set()
for i, l in enumerate(sys.stdin):
    if not i:
        sys.stdout.write(l)
        continue
    ldata = l[:-1].split()
    gid = ldata[0]
    if gid in ids:
        if gid.split('.')[-1].isdigit():
            gid = ".".join(gid.split('.')[:-1])+int(gid.split('.')[-1])+1
        else:
            gid = gid+".1"
    ids.add(gid)
    ndata = [gid,] + [str(int(round(float(x)))) for x in ldata[1:]]
    sys.stdout.write("\t".join(ndata)+"\n")    