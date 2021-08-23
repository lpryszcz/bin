#!/usr/bin/env python2
"""Convert binary ABI sequence format (.ab1) to fastq"""

import sys
from Bio import SeqIO

for fn in sys.argv[1:]:
    r = SeqIO.parse(fn, 'abi').next()
    r.id = r.name
    sys.stdout.write(r.format('fastq'))
    