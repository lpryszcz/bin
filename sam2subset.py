#!/usr/bin/env python
# Select reads from SAM starting withing BED intervals
# and report to seperate output files accordingly to BED name field
# USAGE: sam2subset.py BAM BED

import pysam, sys
bam, bed = sys.argv[1:3]

# open existing bam file
sam = pysam.AlignmentFile(bam)

for l in open(bed):
    ldata = l[:-1].split()
    chrom, s, e, name = ldata[:4]
    s, e = int(s), int(e)
    outfn = "%s.%s.bam"%(bam[:-4], name); print(outfn)
    out = pysam.AlignmentFile(outfn, "wb", template=sam)
    j = k = 0
    for j, a in enumerate(sam.fetch(chrom, s, e), 1):
        if a.pos>s and a.pos<=e:
            out.write(a)
            k += 1
    print(" %s out of %s reads for %s:%s-%s saved"%(k, j, chrom, s, e))
    out.close()
    pysam.index(outfn)
    
    
