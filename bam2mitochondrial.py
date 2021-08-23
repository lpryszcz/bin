#!/usr/bin/env python2
"""
Print contigs that are significantly over-represented in terms of mapped reads.

Usage:
s=c.meta94DNS; samtools idxstats $s.paired.bam | bam2mitochondrial.py -f 5.0 -i ../assemblies/$s.fa -o $s.mitochondrion.fa

"""


import sys
from Bio import SeqIO
from datetime import datetime
from optparse import OptionParser
from numpy    import std, median

def weigthed_std( c2l,c2r ):
  """
  """
  for c in c2l:
    l=c2l
    r=c2r
    cov = round(r*1.0/l,1)

def bam2coverage( handle,fDiff=1.5,iFile=None,outFn=None ):
  """
  """
  fastas = {}
  if iFile:
    fastas = SeqIO.to_dict( SeqIO.parse( open(iFile),'fasta' ) )
  
  outFile = None
  if outFn and fastas:
    outFile = open(outFn,'w' )
  
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
  sd     = std( covs ) #weigthed_std( c2l,c2r ) #
  med    = median( covs )
  
  print "%s reads aligned onto %s bases in %s contigs ( mean = %.2f, sd = %2.f, median = %s)." % ( tReads,tLen,len(c2l),cov,sd,med )
  
  i=0
  diffLen=0
  contigs=[]
  #get contigs starting from highest coverage
  for c in sorted( c2l, key=lambda c: c2r[c]*1.0/c2l[c], reverse=True ):
    l=c2l[c]
    r=c2r[c]
    m=r*1.0/l
    if m>fDiff*cov:# or m<cov/fDiff:
      i+=1
      diffLen+=l
      print "%s\t%s\t%s\t%s\t%.2f" % ( i,c,l,r,m )
      #write contig fasta
      if outFile:
        seq = str( fastas[c].seq )
        outFile.write( '>%s coverage:%.2f; mean:%.2f\n%s\n' % (c,m,cov,seq) )
      
  print "%s bases in %s contigs > %s x mean." % ( diffLen,i,fDiff )#,fDiff )  # OR < %s x mean
  if outFile:
    outFile.close()

def main():
  parser = OptionParser() #allow_interspersed_args=True
  
  #parser.add_option("-v", dest="verbose", action="store_true", default=False,
  #                  help="print status messages to stdout [default: %default]")
  parser.add_option("-f", dest="fDiff", default=1.5, type=float,
                    help="coverage fold change [default: %default]")
  parser.add_option("-i", dest="iFile", default=None, type=str,
                    help="fasta with contig sequences [default: %default]")              
  parser.add_option("-o", dest="oFile", default=None, type=str,
                    help="out file for fasta sequences [default: %default]")
                    
  #read parameters             
  o,args = parser.parse_args()  
  
  handle=sys.stdin
  
  bam2coverage( handle,o.fDiff,o.iFile,o.oFile )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
