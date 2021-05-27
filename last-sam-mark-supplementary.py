#!/usr/bin/env python3
# Mark supplementary alignments in LASTal output.
# this works for maf-convert sam | samtools view -ShT $ref algs
import os, sys, pysam

def store_algs(out, read_algs):
    """Store only reads with mapq lower than max as supplementary algs"""
    if len(read_algs)>1:
        mapq = [a.mapq for a in read_algs]
        maxq = max(mapq)
        # skip cases where more than 1 alg has top mapq
        if len([a for a in read_algs if a.mapq==maxq])>1: maxq+=1
        # mark supplementary
        for a in read_algs:
            if a.mapq<maxq: a.is_supplementary = True
    # report algs
    for a in read_algs: out.write(a)

#samfn = "/home/lpryszcz/cluster/rna_mods/tRNA/_archives/raw/RNA222222.hac.lastal.T0Q1trainD100.maf.sam"
sam = pysam.AlignmentFile(sys.stdin) 
out = pysam.AlignmentFile(sys.stdout,  "wb", template=sam)
pname = None
read_algs = []
for a in sam:
    if pname!=a.qname:
        if read_algs: store_algs(out, read_algs)
        pname = a.qname
        read_algs = []
    # store alg
    read_algs.append(a)
if read_algs: store_algs(out, read_algs)
out.close()
