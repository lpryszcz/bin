#!/usr/bin/env python2
"""
Parse SAM file and output only pairs, in which at least one read from the pair is uniquely map.
It requires XT:A:U and XT:A:R flags in sam alg (ie. BWA output these). 

USAGE: 
samtools view -Shf3 M31.hybrid.sam | python sam2unique.py > M31.hybrid.unique.sam 

M31
HWI-ST227:145:C06RAACXX:7:1101:21063:2636	83	gi|83755449|gb|CP000159.1|	3329954	22	46M	=	3329848	-152	CTTGTCCTGGCAAAGCCAGAACGTTCTCGTGCGACTTGCATGTGTT	GHHCBHJIGIGIJIJIGIGIJGJIGGIIHGIHJJIJIFHHGHFFFA	XT:A:R	NM:i:0	SM:i:0	AM:i:0	X0:i:2	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:46	XA:Z:gi|294342143|emb|FP565814.1|,-3400949,46M,0;
HWI-ST227:145:C06RAACXX:7:1101:21063:2636	163	gi|83755449|gb|CP000159.1|	3329848	23	46M	=	3329954	152	CCCCATGCGGGATCCAGACGACGTGCGGTATTAGCCACGGTTTCCC	FFFFHHHGHJGJJJJJIIIJIJDHIHIIHGFHIIIJFGHHEFFECD	XT:A:U	NM:i:0	SM:i:23	AM:i:0	X0:i:1	X1:i:1	XM:i:0	XO:i:0	XG:i:0	MD:Z:46	XA:Z:gi|294342143|emb|FP565814.1|,+3400843,46M,1;

M8
HWI-ST227:145:C06RAACXX:7:1101:14162:2969	99	gi|294342143|emb|FP565814.1|	3399704	22	46M	=	3399791	133	GCGATTCCTGCTTCATGGAGTCGAGTTGCAGACTCCAATCCGAACT	FFFFHHHGHJJJIJJJJJIIGCGCHFHIJJJBCHGJIJIEIFHGGH	XT:A:R	NM:i:0	SM:i:0	AM:i:0	X0:i:2	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:46	XA:Z:gi|83755449|gb|CP000159.1|,+3328709,46M,0;
HWI-ST227:145:C06RAACXX:7:1101:14162:2969	147	gi|294342143|emb|FP565814.1|	3399791	23	46M	=	3399704	-133	TGCTCGTTGTACCAGCCATTGTAGCACGTGTGCAGCCCTAGGCGTA	8866:?BB9?9JJJJHFIJJJJJJIHJIIIIJJJJJJHHGHHFFFF	XT:A:U	NM:i:0	SM:i:23	AM:i:0	X0:i:1	X1:i:1	XM:i:0	XO:i:0	XG:i:0	MD:Z:46	XA:Z:gi|83755449|gb|CP000159.1|,-3328796,46M,1;

"""
from datetime import datetime
import os, sys

def sam2unique( handle = sys.stdin ):
  """
  """
  i = unique = bothUnique = 0
  uPair = 0
  lines = ''
  for l in handle:
    #write header info
    if l.startswith('@'):
      sys.stdout.write( l )
      continue
      
    #reads
    lines += '%s\n' % '\t'.join( l.split('\t')[:11] )
    i+=1
    if 'XT:A:U' in l:
      uPair += 1
      
    if not i % 2:
      if uPair:
        sys.stdout.write( lines )
        unique += 1
        if uPair == 2:
          bothUnique += 1
          
      uPair = 0
      lines = ''
  
  pairs = i/2    
  sys.stderr.write( 'Processed pairs:\t%s\nUnique pairs:\t%s [%.2f%s]\nBoth unique:\t%s [%.2f%s]\n' % ( pairs,unique,unique*100.0/pairs,'%',bothUnique,bothUnique*100.0/pairs,'%' ) )
    
if __name__=='__main__': 
  T0=datetime.now()
  sam2unique()
  sys.stderr.write( "Elapsed time: %s\n" % ( datetime.now()-T0 ) )
