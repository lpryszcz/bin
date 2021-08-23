#!/usr/bin/env python2
"""Convert to gb/embl to GFF.
Not working
USAGE: toGFF.py in.embl > out.gff
"""

import sys
from BCBio import GFF
from Bio import SeqIO
 
fn = sys.argv[1]
fformat = fn.split('.')[-1]
 
GFF.write(SeqIO.parse(fn, fformat), sys.stdout)
 
