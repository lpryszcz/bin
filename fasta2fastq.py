#!/usr/bin/env python3
# Convert dummy fastq from fasta with default PHRED+33 quality (dummy_qual) of 40 (I). 
# USAGE: cat fasta | fasta2fastq.py [dummy_qual]

import sys

qual = "I"
if len(sys.argv)>2:
    qual = sys.argv[1]

def parse_fasta(handle):
    """Return fasta entries"""
    sid, seqs = "", []
    for l in handle:
        if l.startswith(">"):
            if sid and seqs:
                yield sid, seqs
            sid, seqs = l, []
        else:
            seqs.append(l)
    if sid and seqs:
        yield sid, seqs

for sid, seqs in parse_fasta(sys.stdin):
    sys.stdout.write("@%s%s+\n%s\n"%(sid[1:], "".join(seqs),
                                     qual*(sum(len(x) for x in seqs)-len(seqs))))