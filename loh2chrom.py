#!/usr/bin/env python
"""Generate chromosome characteristics plot for testing for meiosis-associated 
recombination. 

USAGE: loh2chrom.py ../../ref/CANOR.chromosomes.fa.fai Ch_T3_7318.bam.gatk.homo.100bp.cov_flt.bed
"""

import os, sys
import numpy as np
from scipy import stats

ref_fai, bed = sys.argv[1:3]

#load chrom sizes
chrom2size = {}
for l in open(ref_fai):
  chrom, size = l.split()[:2]
  size = int(size)
  chrom2size[chrom] = size

chrom2lohs = {chrom: [] for chrom in chrom2size}  
for l in open(bed):
  chrom, s, e = l.split()[:3]
  s, e = int(s), int(e)
  chrom2lohs[chrom].append(e-s)
  
chrs = [size for chrom, size in sorted(chrom2size.items())]
loh = [np.mean(lohs) for chrom, lohs in sorted(chrom2lohs.items())]
hetero = [100-100.0*sum(lohs)/size for (chrom, size), (chrom, lohs) in zip(sorted(chrom2size.items()), sorted(chrom2lohs.items()))]

#print hetero
#hetero=[18.68,15.37,19.49,26.98,15.01,12.13,9.20,12.58]
#loh=[1980.42,2321.94,2039.60,1214.59,2649.85,3190.51,3765.92,2806.15]
loh=np.array(loh) / 1e3
#chrs=[2936865,2433673,1644167,1591921,1474802,1027593,937312,613068]
chrs=np.array(chrs) / 1e6

print "#sample\tchromosome size\tLOH total\tLOH median\tLOH mean\tLOH stdev"
for (chrom, size), (chrom, lohs) in zip(sorted(chrom2size.items()), sorted(chrom2lohs.items())):
  print "%s\t%s\t%s\t%s\t%s\t%s"%(chrom, size, sum(lohs), np.median(lohs), np.mean(lohs), np.std(lohs))

print "LOH mean size vs Chromosome size"
print " Pearson: r=%s p=%s" % stats.pearsonr(loh, chrs)
print " Spearman: r=%s p=%s" % stats.spearmanr(loh, chrs)
print "Heterozygous % vs Chromosome size"
print " Pearson: r=%s p=%s" % stats.pearsonr(hetero, chrs)
print " Spearman: r=%s p=%s" % stats.spearmanr(hetero, chrs)

import matplotlib.pyplot as plt

fig = plt.figure()
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
p1, = ax1.plot(chrs, hetero, "bo", label="Heterozygous")
ax1.set_ylabel("Heterozygous [%]")
ax2.set_xlabel("Chromosome size [Mb]")
p2, = ax2.plot(chrs, loh, "ro", label="LOH mean size")
ax2.set_ylabel("Mean size [kb]")
plt.show()



