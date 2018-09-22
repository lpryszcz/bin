#!/usr/bin/env python
"""print common SNP info

USAGE: snp2common.py SNPs1 SNPs2 [...] SNPsN
"""

import sys
import numpy as np

fnames = sys.argv[1:]

import scipy
import pylab
import scipy.cluster.hierarchy as sch


#load
sys.stderr.write("Parsing...\n")
samples = []
s2snps = {}
for i, fn in enumerate(fnames, 1):
	s = fn.split(".bam")[0]
	sys.stderr.write(" %s / %s %s %s                   \r"%(i, len(fnames), fn, s))
	samples.append(s)
	if s not in s2snps:
		s2snps[s] = set()
	for l in open(fn):
		if l.startswith(('#',"scaffold6|")):
			continue
		ldata = l.split('\t')
		#bam2snp output
		if ':' in ldata[0]:
			coord, refc, refb, refFreq, cov, base, freq = ldata[:7]
		#vcf
		else:
			chrom, pos, name, refb, base, qual, filt = ldata[:7]
			coord = "%s:%s"%(chrom, pos)
		s2snps[s].add((coord, refb, base))

#report
sys.stderr.write("Computing distances...\n")
a = np.zeros((len(s2snps), len(s2snps)))
for i in range(len(s2snps)):
	s1 = samples[i]
	for j in range(i, len(s2snps)):
		s2 = samples[j]
		#a[i][j] = a[j][i] = len(s2snps[s1].intersection(s2snps[s2]))
		a[i][j] = a[j][i] = len(s2snps[s1].symmetric_difference(s2snps[s2]))
		
np.savetxt(sys.stdout, a, delimiter="\t", header="\t".join(samples[:len(s2snps)]), fmt='%i')

sys.stderr.write("Plotting...\n")
D=a
# Compute and plot first dendrogram.
fig = pylab.figure(figsize=(8,8))
ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
Y = sch.linkage(D, method='centroid')
Z1 = sch.dendrogram(Y, orientation='left')
ax1.set_xticks([])
#ax1.set_yticks([])

# Compute and plot second dendrogram.
ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
Y = sch.linkage(D, method='single')
Z2 = sch.dendrogram(Y)
#ax2.set_xticks([])
ax2.set_yticks([])

fig.savefig('dendrogram.svg')

# Plot distance matrix.
axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
idx1 = Z1['leaves']
idx2 = Z2['leaves']
D = D[idx1,:]
D = D[:,idx2]
im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
axmatrix.set_xticks([])
axmatrix.set_yticks([])

# Plot colorbar.
axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
pylab.colorbar(im, cax=axcolor)
fig.show()
fig.savefig('dendrogram.full.svg')	

