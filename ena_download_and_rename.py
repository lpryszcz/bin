#!/usr/bin/env python3
# ./download.py ../PRJNA300878.txt [name_field www_field]

import os, sys

fname = sys.argv[1]
name_field, www_field = 9, 10
if len(sys.argv)>3:
    name_field, www_field = map(int, sys.argv[2:4])

for l in open(fname):
    if l.startswith('study_accession'): continue
    ldata = l[:-1].split('\t'); print(ldata)
    name, wwws = ldata[name_field], ldata[www_field].split(';')
    for pi, www in enumerate(wwws, 1):
        paired = ""
        if len(wwws)>1:
            paired="_%s"%pi
        cmds = "wget -nc %s && ln -s %s %s%s.fq.gz"%(www, os.path.basename(www), name, paired); print(cmds)
        os.system(cmds)
    #break
