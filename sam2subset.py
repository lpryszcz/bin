#!/usr/bin/env python2
# Select reads from SAM starting within BED intervals +/- n [5]
# and report to seperate output files accordingly to BED name field
# USAGE: sam2subset.py BAM BED [n]

import pysam, sys
bam, bed = sys.argv[1:3]
n = 5

# open existing bam file
sam = pysam.AlignmentFile(bam)

for l in open(bed):
    ldata = l[:-1].split()
    chrom, s, e, name = ldata[:4]
    # get minlen
    if len(ldata)>4 and ldata[4].isdigit():
        minlen = int(ldata[4])
    else:
        minlen = 0
    s, e = int(s), int(e)
    s = s - n if s>=n else 0
    e = e + n
    outfn = "%s.%s.bam"%(bam[:-4], name); print(outfn)
    out = pysam.AlignmentFile(outfn, "wb", template=sam)
    j = k = 0
    for j, a in enumerate(sam.fetch(chrom, s, e), 1):
        if a.pos>=s and a.pos<=e:
            out.write(a)
            k += 1
    # store those longer than minlen as genomic
    if minlen:
        for a in sam.fetch(chrom, e):
            if a.pos>e and a.alen>minlen:
                out.write(a)
                k += 1
    print(" %s out of %s reads for %s:%s-%s saved"%(k, j, chrom, s, e))
    out.close()
    pysam.index(outfn)
    
    
