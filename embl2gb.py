#!/usr/bin/env python2
"""
Download genomes from NCBI for given taxid.
"""

import sys
from datetime import datetime
from Bio      import SeqIO

def split_seq( seq,length=60 ):
  """Return list of sequence fragnments having given length.
  """
  return [seq[i:i + length] for i in range(0, len(seq), length)]

t0=datetime.now()

records = SeqIO.parse( sys.stdin,'embl' )

SeqIO.write( records,sys.stdout,'gb' )

dt=datetime.now()-t0
sys.stderr.write( "#Time elapsed: %s\n" % dt )
