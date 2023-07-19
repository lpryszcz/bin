#!/usr/bin/env python
# Report N50 from contig lenghts
# n50 of fasta file
# > samtools faidx FASTA; cat FASTA.fai | len2n50.py
# n50 of mapped reads (excluding secondary, supplementary & QC fail reads)
# > samtools view -F 3844 BAM | awk '{print length($10)}' | len2n50.py
# n50 of unmapped reads
# > samtools view -f 4 BAM | awk '{print length($10)}' | len2n50.py

import sys

def get_N50(lengths, genomeFrac=0.5):
    """Return N50 from contig lengths"""
    genomeSize = sum(lengths)
    # parse contigs by descending size
    size = totsize = 0
    for i, size in enumerate(sorted(lengths, reverse=True), 1):
        totsize += size
        if totsize >= genomeFrac*genomeSize:
            break
    return size

lengths = [int(l.strip()) for l in sys.stdin]
n50 = get_N50(lengths)
print(n50)

