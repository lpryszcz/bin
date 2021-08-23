#!/usr/bin/env python2

"""
Usage:
samtools idxstats CBS1954.paired.bam | bam2coverage.py 

"""


import sys
#from Bio import SeqIO
from datetime import datetime
from numpy import std,median

def bam2coverage():
  """
  """
  fDiff=1.5
  if len( sys.argv )>1: fDiff=float( sys.argv[1] )
  
  handle=sys.stdin
  c2l,c2r,c2m={},{},{}
  covs=[]
  for line in handle:
    c,l,r,unmapped=line.split('\t')
    l,r=int(l),int(r)
    if l<1000: continue
    c2l[c]=l
    c2r[c]=r
    covs.append( round(r*1.0/l,1) ) # 123.333 -> 123.3
  
  tLen   = sum( c2l.itervalues() )
  tReads = sum( c2r.itervalues() )
  cov    = tReads*1.0/tLen
  sd     = std( covs )
  med    = median( covs )
  
  print "%s reads aligned onto %s bases in %s contigs ( mean = %.2f, sd = %2.f, median = %s)." % ( tReads,tLen,len(c2l),cov,sd,med )
  
  i=0
  diffLen=0
  for c in c2l:
    l=c2l[c]
    r=c2r[c]
    m=r*1.0/l
    if m>fDiff*cov:# or m<cov/fDiff:
      i+=1
      diffLen+=l
      #print "%s\t%s\t%s\t%s\t%.2f" % ( i,c,l,r,m )
  print "%s bases in %s contigs > %s x mean." % ( diffLen,i,fDiff )#,fDiff )  # OR < %s x mean
  

if __name__=='__main__': 
  t0=datetime.now()
  bam2coverage()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
