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
        # add children
        t.add_child(i2n[idx1], dist=dist)
        t.add_child(i2n[idx2], dist=dist)
        i2n[n+idx] = t
        idx += 1
    return t 

def get_name(contig):
    return contig.split('|')[-1].split('.')[0]
            
def array2clusters(infile):
    # load matrix
    tr.logger("Loading matrix from %s ..."%infile)
    d, bin_chr, bin_position = tr.load_data_txt(infile, remove_nans=True, chrs=[], retain=1, remove_shorter=0)
    genomeSize = np.diff(bin_position, axis=1).sum()
    contig2size = {get_name(c): 0 for c in np.unique(bin_chr)}
    for c, (s, e) in zip(bin_chr, bin_position):
        contig2size[get_name(c)] += e-s
    print " loaded %s contigs summing %s bp"%(d.shape[0], genomeSize)
        
    # cluster
    tr.logger("Calculating linkage matrix...")
    method="ward"
    transform = lambda x: np.log(np.max(x+1))-np.log(x+1)
    Z = sch.linkage(transform(d)[np.triu_indices(d.shape[0], 1)], method=method); del d
    
    # get ete Tree
    tr.logger("Populating tree...")
    names = ["%s.%s"%(get_name(c), s/1000) for c, (s, e) in zip(bin_chr, bin_position)]
    t = array2tree(Z, names); del Z
    
    # generate clusters
    tr.logger("Generating subtrees...")
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
    
    tr.logger("Assigning contigs to %s clusters..."%len(subtrees))
    total = correct = 0
    contig2cluster = {get_name(c): Counter() for c in np.unique(bin_chr)}
    for i, subtree in enumerate(subtrees, 1): 
        c = Counter(get_name(_n.name) for _n in subtree)
        #print " %s %s %s"%(i, len(subtree), c.most_common(3))
        total += len(subtree)
        correct += c.most_common(1)[0][1]
        # poplate contig2clustre
        for k, v in c.iteritems():
            if not k: continue
            contig2cluster[k][i] += v
    print " %s / %s [%.2f%s]"%(correct, total, 100.*correct/total, '%')

    tr.logger("Weak assignments...")
    clusters = [[] for i in range(len(subtree))]
    withoutCluster, weakCluster = [], []
    for c, counter in contig2cluster.iteritems():
        if not counter:
            withoutCluster.append(c)
            continue
        # get major cluster
        clusteri, count = counter.most_common(1)[0]
        mfrac = 1. * count / sum(counter.itervalues())
        clusters[clusteri].append(c)
        if mfrac<.66:
            #print " %s %s: %s" %(c, sum(counter.itervalues()), str(counter.most_common(3)))
            weakCluster.append(c)
    print " %s bp in %s contigs without assignment."%(sum(contig2size[c] for c in withoutCluster), len(withoutCluster))
    print " %s bp in %s contigs having weak assignment."%(sum(contig2size[c] for c in weakCluster), len(weakCluster))
 
    clusters = filter(lambda x: x, clusters)
    outfile = infile+".clusters.tab"
    tr.logger("Reporting clusters to: %s ..."%outfile)
    totsize = 0
    with open(outfile, "w") as out:
        for i, cluster in enumerate(clusters, 1):
            clSize = sum(contig2size[c] for c in cluster)
            totsize += clSize
            print " %s %s bp in %s contigs" % (i, clSize, len(cluster))
            out.write("\t".join(cluster)+"\n")
    print "%3s bp in %s clusters!"%(totsize, len(clusters))

def main():
    infile=sys.argv[1] 
    array2clusters(infile)
    
if __name__=="__main__":
    t0 = datetime.now()
    main()
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
