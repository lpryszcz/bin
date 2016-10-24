#!/usr/bin/env python

import os, sys
import scipy.cluster.hierarchy as sch
import ete3
import numpy as np
sys.path.append("/home/lpryszcz/src/dna-triangulation/")
import triangulation as tr
from collections import Counter
from datetime import datetime

def array2tree(Z, names=[]):
    """Return tree representation for array"""
    n = Z.shape[0]+1
    if not names:
        names = map(str, range(n))
        
    i2n = {}
    t = ete3.Tree()
    idx = 0
    for i, (idx1, idx2, dist, sample_count) in enumerate(Z, 1):
        idx1, idx2 = int(idx1), int(idx2)
        if idx1<n:
            i2n[idx1] = ete3.Tree(name=names[idx1], dist=0)
        if idx2<n:
            i2n[idx2] = ete3.Tree(name=names[idx2], dist=0)
        t = ete3.Tree(dist=0)
        # normalise distance
        dist -= max(i2n[idx1].get_farthest_leaf()[1], i2n[idx2].get_farthest_leaf()[1])
        # print i, idx1, idx2, dist, ddist, rthest_leaf()[1]
        t.add_child(i2n[idx1], dist=dist)#, support=sample_count)
        t.add_child(i2n[idx2], dist=dist)#, support=sample_count)
        i2n[n+idx] = t
        idx += 1

    return t #, i2n

def array2clusters(infile):
    # load matrix
    d, bin_chr, bin_position = tr.load_data_txt(infile, remove_nans=True, chrs=[], retain=1, remove_shorter=0)
    genomeSize = np.diff(bin_position, axis=1).sum()
    print "loaded %s contigs summing %s bp"%(d.shape[0], genomeSize)
    # cluster
    method="ward"
    transform = lambda x: np.log(np.max(x+1))-np.log(x+1)
    Z = sch.linkage(transform(d)[np.triu_indices(d.shape[0], 1)], method=method)
    # get ete Tree
    names = ["%s.%s"%(c.split('|')[-1], s/1000) for c, (s, e) in zip(bin_chr, bin_position)]
    t = array2tree(Z, names)#; print t.describe()#; t.show() #'''
    # generate clusters
    subtrees=[]
    while len(t)>2:
        n, dist = t.get_farthest_leaf()
        dists = [a.dist for a in n.get_ancestors()]
        # get ancestor with the longest branch length
        ai = dists.index(max(dists))
        a = n.get_ancestors()[ai]
        if n.name:
            c = Counter(_n.name.split('.')[0] for _n in a)
            subtrees.append(a)
        p = n.get_ancestors()[ai+1]
        p.remove_child(a)
    
    print "Formed %s clusters."%len(subtrees)
    total = correct = 0
    for i, a in enumerate(subtrees, 1):
        c = Counter(_n.name.split('.')[0] for _n in a)
        print " %s %s %s"%(i, len(a), c.most_common(10))
        total += sum(c.values())
        correct += c.most_common(1)[0][1]
    print "%s / %s [%.2f%s]"%(correct, total, 100.*correct/total, '%')

def main():
    infile=sys.argv[1] 
    array2clusters(infile)
    
if __name__=="__main__":
    t0 = datetime.now()
    main()
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
