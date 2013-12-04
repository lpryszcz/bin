#!/usr/bin/env python
"""
"""
import os,sys
import numpy as np
from datetime import datetime
from math import fabs

def main( args=[] ):
  maxInsert=999999
  if args: maxInsert=int(args[0])
  
  handle=sys.stdin
  insertsSizes=[]
  for line in handle:
    if line.startswith("@"):
      continue
    QNAME,FLAG,RNAME,POS,MAPQ,CIAGR,MRNM,MPOS,TLEN,SEQ,QUAL=line.split('\t')[:11]
    TLEN=fabs( int(TLEN) )#; print TLEN
    if TLEN and TLEN<maxInsert: insertsSizes.append( TLEN )
    
  mu = np.mean( insertsSizes )
  sd = np.std( insertsSizes )
  median = np.median( insertsSizes )
  #print "#mean\tsd\tmedian"
  print "%s\t%s\t%s" % (mu,sd,median)
  
if __name__=='__main__': 
  t0=datetime.now()
  main( sys.argv[1:] )
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
