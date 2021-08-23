#!/usr/bin/env python2
"""
Parse SAM file and output only pairs uniquely aligned.

USAGE: 
samtools view -St yeast_chromosomes.fa.fai 409.sam -f3 | samgem2unique.py > 409.unique.sam 


"""
from datetime import datetime
import os, sys

def report(refs,lines):
  """ """
  if lines and len(refs)==2 and refs[0]==refs[1]:
    sys.stdout.write( lines )
    return 1
  return 0

def sam2unique( handle = sys.stdin ):
  """Add first/second in pair checking!"""
  i = k = unique = bothUnique = 0
  pName = lines = ''
  refs=[]
  for l in handle:
    #write header info
    if l.startswith('@'):
      sys.stdout.write( l )
      continue
      
    #name,flag,contig,pos,mapq,cigar,paired,pairStart,isize,seq,qual
    data   = l.split('\t')
    
    name,flag,ref = data[:3]
    if name != pName:
      i+=1
      unique += report(refs,lines)
      #reset variables
      lines = l
      refs = [ ref ]
      pName = name
    else:
      #reads
      refs.append( ref )
      lines += l
      
  #report
  unique += report(refs,lines)
  
  sys.stderr.write( 'Processed pairs:\t%s\nUnique pairs:\t%s [%.2f%s]\n' % ( i,unique,unique*100.0/i,'%' ) ) #,bothUnique,bothUnique*100.0/pairs,'%'
    
if __name__=='__main__': 
  T0=datetime.now()
  sam2unique()
  sys.stderr.write( "Elapsed time: %s\n" % ( datetime.now()-T0 ) )
