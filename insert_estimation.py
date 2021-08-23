#!/usr/bin/env python2
"""
Estimate insert size.

USAGE:
insert_estimation.py ref fastq1 fastq2

optionally - build crappy assembly [not yet]
make alignment
get insert distribution

Author:
lpryszcz@crg.es

Required:
samtools
bowtie
SOAPdenovo (optionally)
"""

import os, sys
import numpy as np
from commands import getoutput
from datetime import datetime
from math import fabs

def insertSize( handle,maxInsert=100000 ):
  """
  """
  insertsSizes=[]
  for line in handle:
    QNAME,FLAG,RNAME,POS,MAPQ,CIAGR,MRNM,MPOS,TLEN,SEQ,QUAL=line.split('\t')[:11]
    TLEN=fabs( int(TLEN) )#; print TLEN
    if TLEN and TLEN<maxInsert: 
      insertsSizes.append( TLEN )
    
  mu = np.mean( insertsSizes )
  sd = np.std( insertsSizes )
  median = np.median( insertsSizes )
  sys.stderr.write( "#mean\tsd\tmedian\n" )
  sys.stderr.write( "%.2f\t%.2f\t%.2f\n" % (mu,sd,median) )
  print "%.2f\t%.2f\t%.2f" % (mu,sd,median)

def main():
  if len(sys.argv)<3:
    sys.exit( "insert_estimation.py ref.fa fastq1 fastq2 [#reads_limit]" )
  ref,q1,q2 = sys.argv[1:4]
  l = 10**5
  if len( sys.argv ) > 4:
    l = int( sys.argv[4] )
  
  #sys.stderr.write( "Selecting %s reads...\n" % ( l, ) )
  #for q in ( q1,q2 ):
  # getoutput( "head -n%s %s > %s.tmp" % (l*4,q,q,) )
  
  if not os.path.isfile( ref+'.1.ebwt' ):
    sys.stderr.write( "Generating bowtie index...\n" )
    getoutput( "bowtie-build %s %s" % (ref,ref) )
  
  sys.stderr.write( "Aligning...\n" )
  cmd1="bowtie --sam --phred33-quals --maxins 50000 -p 1 -q -u %s --fr %s -1 %s -2 %s | samtools view -f3 -S - > %s.sam" % ( l,ref,q1,q2,q1 )
  handle=getoutput( cmd1 )
  
  sys.stderr.write( "Calculating insert disrtibution...\n" )
  handle = open( "%s.sam" % q1 )
  insertSize( handle )
  
  sys.stderr.write( "You can now remove bowtie-index files: rm *.ebwt %s.sam\n" % ( q1, ) )
  
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )

