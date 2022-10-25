#!/usr/bin/env python
# Report N50 from contig lenghts

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

