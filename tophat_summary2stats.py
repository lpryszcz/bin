#!/usr/bin/env python2
#Report alignement stats from tophat2 summary
#USAGE: cat align_summary.txt | python tophat_summary2stats.py

import sys

def get_stats(fn):
    """Return stats"""
    reads = algs = multi_algs = corcordant = uniq_corcordant = 0
    for l in open(fn):
        ldata = l[:-1].split()
        if not ldata:
            continue
        #          Input     :   6660317
        if ldata[0] ==   'Input':
            reads += int(ldata[2])
        #           Mapped   :   5864364 (88.0% of input)            
        elif ldata[0] == 'Mapped':
            algs += int(ldata[2])
        #Aligned pairs:   5305465
        elif ldata[0] == 'Aligned': 
            corcordant += int(ldata[2])
        elif ldata[0] == 'of':
            #     of these:     84764 ( 1.6%) have multiple alignments        
            if corcordant:
                uniq_corcordant = corcordant - int(ldata[2])
            #     of these:     96376 ( 1.6%) have multiple alignments (1 have >20)                
            else:
                multi_algs += int(ldata[2])
        elif ldata[0].isdigit():
            corcordant -= int(ldata[0])                
    return reads, algs, algs-multi_algs, corcordant

print "#sample\treads\taligned\t[%]\taligned uniqly\t[%]\tcorcordant pairs\t[%]"
outline = "%s\t%s\t%s\t%.2f%s\t%s\t%.2f%s\t%s\t%.2f%s"
for fn in sys.argv[1:]:
    reads, algs, uniq, corcordant = get_stats(fn)
    print outline % (fn, reads, algs, algs*100.0/reads, '%', \
                     uniq, uniq*100.0/reads, '%', \
                     corcordant, corcordant*200.0/reads, '%')
