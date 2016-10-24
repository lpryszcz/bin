#!/usr/bin/env python

import os, sys
import scipy.cluster.hierarchy as sch
import ete3
import numpy as np
sys.path.append("/home/lpryszcz/src/dna-triangulation/")
import triangulation as tr
from collections import Counter


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

def main(infile="../_archives/snap/SRR2626163.100k.npz")
    d, bin_chr, bin_position = tr.load_data_txt(infile, remove_nans=True, chrs=[], retain=1, remove_shorter=0)
    genomeSize = np.diff(bin_position, axis=1).sum()
    print "loaded %s contigs summing %s bp"%(d.shape[0], genomeSize)
    
    method="ward"
    transform = lambda x: np.log(np.max(x+1))-np.log(x+1)
    Z = sch.linkage(transform(d)[np.triu_indices(d.shape[0], 1)], method=method)
    names = ["%s.%s"%(c.split('|')[1], s/1000) for c, (s, e) in zip(bin_chr, bin_position)]
    '''
    while tr._check_hierarchy_uses_cluster_more_than_once(Z):
        wrong = tr._check_hierarchy_uses_cluster_more_than_once(Z); print d.shape, Z.shape, wrong
        d = np.delete(np.delete(d, wrong, 1), wrong, 0)
        bin_chr = np.delete(bin_chr, wrong, 0)
        bin_position= np.delete(bin_position, wrong, 0)
        sys.stderr.write("[WARNING] Ward failed. Dropping %s problematic clusters...\n"%len(wrong))
        #Z = np.delete(Z, wrong, 0)
        Z = sch.linkage(transform(d)[np.triu_indices(d.shape[0], 1)], method=method)

    dn = sch.dendrogram(Z)
    tree = sch.to_tree(Z)
    nw = tr.getNewick(tree, "", tree.dist, names)

    t = ete3.Tree(nw); print t.describe()#; t.show()
    ''' 
    t = array2tree(Z, names); print t.describe()#; t.show() #'''
    subtrees=[]
    while len(t)>2:
        n, dist = t.get_farthest_leaf()
        dists = [a.dist for a in n.get_ancestors()]
        #print len(t), n.name, dists
        #idist = [] #i for i in range(len(dists)-1) if dists[i]>2*dists[i+1]]
        idist = []# [i for i in range(1, len(dists)) if dists[i]>2*dists[i-1]]
        if idist:
            ai = idist[0]
        else:
            ai = dists.index(max(dists))
        a = n.get_ancestors()[ai] 
        if n.name:
            c = Counter(_n.name.split('.')[0] for _n in a)
            subtrees.append(a)
            #print " ", a.dist, c.most_common(10)
        p = n.get_ancestors()[ai+1]
        p.remove_child(a)#; print t.describe()

    print len(subtrees)
    total = correct = 0
    for i, a in enumerate(subtrees, 1):
        c = Counter(_n.name.split('.')[0] for _n in a)
        print i, len(a), c.most_common(10)
        total += sum(c.values())
        correct += c.most_common(1)[0][1]
    print "%s / %s [%.2f%s]"%(correct, total, 100.*correct/total, '%')
    
if __name__=="__main__":
    infile=sys.argv[1] 
    t0 = datetime.now()
    main(infile)
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
