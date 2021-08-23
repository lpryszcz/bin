#!/usr/bin/env python2
# USAGE: ./bed2cnv_plot.py cnv.100.bed Title


import sys
# http://matplotlib.org/examples/api/two_scales.html
import numpy as np
import matplotlib.pyplot as plt

def load_bed(fnames, minSize=0):
  """Return chrs and sizes of cnvs within each chromosome"""
  chrom2cnvs = {}
  for fn in fnames:
    for l in open(fn):
      if l.startswith('#'): continue
      chrom, s, e = l[:-1].split()[:3]
      s, e = int(s), int(e)
      size = e-s
      if size<minSize: continue
      if chrom not in chrom2cnvs:
        chrom2cnvs[chrom] = []
      chrom2cnvs[chrom].append(size)
  # get sorted chromosomes and cnv sizes
  chrs = sorted(filter(lambda x: x.isdigit(), chrom2cnvs), key=lambda x: int(x)) # .startswith('DS')
  cnvs = [chrom2cnvs[chrom] for chrom in chrs]
  # add the rest
  chrs.append("*")
  cnvs.append([])
  for chrom in filter(lambda x: x not in chrs, chrom2cnvs):
    cnvs[-1]+= chrom2cnvs[chrom]
  return chrs, cnvs

fnames = sys.argv[1:-1]
title = sys.argv[-1]

chrs, cnvs = load_bed(fnames)

width, spacer = 0.4, 0.1
color1, color2 = 'b', 'r'

ind = np.arange(len(chrs))#; print chrs, ind

fig, ax1 = plt.subplots()
ax1.bar(ind+spacer, map(len, cnvs), width, color=color1)
ax1.set_xlabel('chromosome')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('number of CNVs', color=color1)
for tl in ax1.get_yticklabels(): tl.set_color(color1)

ax2 = ax1.twinx()
ax2.bar(ind+width+spacer, np.array(map(sum, cnvs))/1e6, width, color=color2)
ax2.set_ylabel('cummulative size of CNVs [Mb]', color=color2)
for tl in ax2.get_yticklabels(): tl.set_color(color2)

ax1.set_xticks(ind+width+spacer)
ax1.set_xticklabels(chrs)

ax1.set_title(title)

ext = ".png"
if len(fnames)==1:
  fname = fnames[0]+ext
else:
  fname = "all%s"%ext
print "Saving figure as %s"%fname
plt.savefig(fname)


