#!/usr/bin/env python
# Produce venn diagrams. So far only reporting overlapping entries. 
# USAGE: venn.py file1 file2 [file3 ... fileN]

import os, sys

if len(sys.argv)<3:
    sys.exit("Provide 2 or more files!\n ie. venn.py file1 file2 [file3 ... fileN]\n")
    
fnames = sys.argv[1:]; print "%s input files: %s"%(len(fnames)," ".join(fnames))

def get_id(l):
    """Return id set. EnsemblIDs are stored as int."""
    return l.split('\t')[0].strip()

starter, linker = "venn-", "-"
venns, names = [], []

# load ids and get overlaps
for i, fn in enumerate(fnames):
    ids = set(get_id(l) for l in open(fn) if not l.startswith('#'))
    print i+1, fn, len(ids)#, list(ids)[:10]
    nvenns, nnames = [], []
    for j, s in enumerate(venns):
        nvenns.append(ids.intersection(s))
        nnames.append(names[j]+linker+fn)
    # insert fn and ids
    names.insert(i, starter+fn)
    venns.insert(i, ids)
    venns += nvenns
    names += nnames
    
# remove overlaps
for i, n in enumerate(names[len(fnames):], len(fnames)):
    # get fnames idx that have to be excluded from given group
    notin = [fnames.index(fn) for fn in set(fnames).difference(n.split(linker))]#; print n, notin
    if notin:
        venns[i] = venns[i].difference(set().union(*[venns[j] for j in notin]))
    
# get uniq per group
for i in range(len(fnames)):
    venns[i] = venns[i].difference(set().union(*venns[:i]+venns[i+1:]))

#print len(venns)
if not os.path.isdir('venn'): os.makedirs('venn')

print "\n%s groups: "%len(venns)
for i, (n, s) in enumerate(zip(names, venns), 1):
    print '', i, len(s), n
    # report combined lines
    header = ''
    id2txt = {x: [] for x in s}
    for j, fn in enumerate(n[len(starter):].split(linker)):
        for k, l in enumerate(open(fn)):
            if l.startswith('#'):
                header += l[:-1]
                continue
            if get_id(l) in id2txt:
                id2txt[get_id(l)].append(l[:-1])
    # 
    with open(n+".txt", 'w') as out:
        out.write(header+'\n'+'\n'.join('\t'.join(x) for x in id2txt.itervalues())+'\n')


