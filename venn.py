#!/usr/bin/env python
# Calculate venn diagrams

import sys

fnames = sys.argv[1:]; print "%s input files: %s"%(len(fnames)," ".join(fnames))

def get_ids(fn):
    """Return id set. EnsemblIDs are stored as int."""
    ids = set()
    for l in open(fn):
        if l.startswith('#'):
            continue
        eid = int(l.split('\t')[0].strip()[7:])
        ids.add(eid)
    return ids

venns, names = [], []

# load ids and get overlaps
for i, fn in enumerate(fnames):
    #ids = get_ids(fn)
    ids = set(l.split('\t')[0].strip() for l in open(fn) if not l.startswith('#'))
    print i+1, fn, len(ids)#, list(ids)[:10]
    nvenns = []
    for j, s in enumerate(venns):
        nvenns.append(ids.intersection(s))
        names.append(names[j]+" - "+fn)
    # insert fn and ids
    names.insert(i, fn)
    venns.insert(i, ids)
    venns += nvenns
    
# get uniq per group
for i in range(len(fnames)):
    venns[i] = venns[i].difference(set().union(*venns[:i]+venns[i+1:]))

#print len(venns)
print "\n%s groups: "%len(venns)
for i, (n, s) in enumerate(zip(names, venns), 1):
    print i, n, len(s)

