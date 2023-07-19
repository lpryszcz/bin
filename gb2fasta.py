#!/usr/bin/env python3
"""
Genbank to FASTA converter.

USAGE:
cat infile | gb2fasta.py > outfile

Author:
l.p.pryszcz+git@gmail.com
"""

import sys
from datetime import datetime
from Bio      import SeqIO

def split_seq( seq,length=60 ):
  """Return list of sequence fragments having given length."""
  return [seq[i:i + length] for i in range(0, len(seq), length)]

t0=datetime.now()

for gb in SeqIO.parse( sys.stdin,'gb' ):
  seq = '\n'.join( split_seq(str(gb.seq)) )
  fasta = ">%s\n%s\n" % ( gb.name,seq )
  sys.stdout.write( fasta )

dt=datetime.now()-t0
sys.stderr.write( "#Time elapsed: %s\n" % dt )
