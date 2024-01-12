#!/usr/bin/env python3
# Select at most N longest reads from transcripts in BAM
# and save them in OUTBAM
# USAGE: bam2select.py BAM OUTBAM N

import pysam, sys

mapq = 15
bam, outbam, n = sys.argv[1:4]
n = int(n)

sam = pysam.AlignmentFile(bam)
out = pysam.AlignmentFile(outbam, "wb", reference_names=sam.references, reference_lengths=sam.lengths)
log = sys.stderr

refs = sam.references
for i, ref in enumerate(refs, 1):
    if not i%10: log.write(" %s / %s %s  \r"%(i, len(refs), ref))
    _n = 0
    for a in sam.fetch(ref):
        # check quality & strand
        if a.mapq<mapq or a.is_secondary or a.is_reverse: continue
        out.write(a)
        _n += 1
        if _n==n: break
                          
out.close()
pysam.index(outbam)
