#!/usr/bin/env python2
"""
Parse SAM file and output only pairs with at least one read aligned.

Compatible with bowtie/bwa output - one entry per read. 
SAM file has to be sorted by read name.

USAGE: 
samtools view -St yeast_chromosomes.fa.fai 409.sam -f3 | sam2aligned.py > 409.aligned.sam 


"""
from datetime import datetime
import os, sys

def int2bin( n, count=12 ):
  """returns the binary of integer n, using count number of digits
  @ http://www.daniweb.com/software-development/python/code/216539
  """
  return "".join([str((n >> y) & 1) for y in range(count-1, -1, -1)])

def sam2unique( handle = sys.stdin ):
  """
  """
  i = k = aligned = 0
  pName = lines = ''
  refs=0
  for l in handle:
    #write header info
    if l.startswith('@'):
      sys.stdout.write( l )
      continue
      
    #name,flag,contig,pos,mapq,cigar,paired,pairStart,isize,seq,qual
    name,flag,ref = l.split('\t')[:3]
    if name != pName:
      i+=1
      if lines and refs:
        sys.stdout.write( lines )
        aligned += 1
        refs = 0
      lines = l
      if ref != "*":
        refs += 1
      pName = name
    else:
      #reads
      if ref != "*":
        refs += 1
      lines += l
      
  if lines and refs:
    aligned += 1
    sys.stdout.write( lines )
  
  sys.stderr.write( 'Processed pairs:\t%s\nAligned pairs:\t%s [%.2f%s]\n' % ( i,aligned,aligned*100.0/i,'%' ) ) #,bothUnique,bothUnique*100.0/pairs,'%'
    
if __name__=='__main__': 
  T0=datetime.now()
  sam2unique()
  sys.stderr.write( "Elapsed time: %s\n" % ( datetime.now()-T0 ) )
