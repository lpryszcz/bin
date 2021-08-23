#!/usr/bin/env python2
"""
Identify SNP sites in mpileup out from BAM alignments. 
In addition calculates overall and each contig coverage. 
 
At first run: /users/tg/lpryszcz/cluster/assembly/mpileup.sh (s=MCO456; samtools mpileup -f gem/$s.fa gem/${s}.bam > gem/${s}.mpileup)
Execute: 
 python ~/workspace/assembly/src/heterozygosity.py -i gem/CBS6318.mpileup [-f minFreq -d minDepth -o out_base_name]
or use with piping:
 s=CBS1954; samtools mpileup -f gem/$s.fa gem/$s.bam | python ~/workspace/assembly/src/heterozygosity.py -o gem/$s -f 0.2 -d 10
"""

import os, sys
from optparse import OptionParser
from datetime import datetime

def _get_base_counts_old( ref,alts,alphabet ):
  """Return a list of A,C,G and T frequencies in given position of the alignment.
  """
  alts=alts.upper()
  #remove indels info
  for symbol in ('-','+'):
    while symbol in alts:
      i=alts.index(symbol)
      baseNo=int(alts[i+1])
      alts=alts[:i]+alts[i+baseNo+2:]
  #count base occurencies
  base_counts=[]
  for base in alphabet:
    if base!=ref: base_counts.append( alts.count(base) ) #base_counts.append( alts.count('.')+alts.count(',') )      
    else:         base_counts.append( 0 )
  return base_counts

def _remove_indels( alts ):
  """
  Remove indels from mpileup.
  .$....,,,,....,.,,..,,.,.,,,,,,,....,.,...,.,.,....,,,........,.A.,...,,......^0.^+.^$.^0.^8.^F.^].^],
  ........,.-25ATCTGGTGGTTGGGATGTTGCCGCT..
  """
  #remove indels info
  for symbol in ('-','+'):
    baseNo = 0
    while symbol in alts:
      i=alts.index(symbol)
      
      j = 1
      digits=[]
      while alts[i+j].isdigit():
        digits.append( alts[i+j] )
        j += 1
      
      if digits:
        baseNo=int( ''.join(digits) )
        
      alts=alts[:i]+alts[i+baseNo+len(digits)+1:] #......+1A..,
      
  return alts

def _get_base_counts( ref,alts,alphabet,illegal=('S','*','[',']','>','<') ):
  """Return a list of A,C,G and T frequencies in given position of the alignment.
  ******CCCCCcccccCccCCCCcccCcCccCcc^SC^]c
  * - deletion
  +INT / -INT - insertion
  S - ??
  [ or ] - ??
  """
  #remove deletions
  dels = alts.count('*') 
  alts = alts.upper()
  #remove insertions
  alts = _remove_indels( alts )
      
  #count base occurencies
  base_counts=[]
  for base in alphabet:
    if base!=ref: base_counts.append( alts.count(base) ) #      
    else:         base_counts.append( alts.count('.')+alts.count(',') )
  return base_counts,dels

def get_alt_allele( base_ref,cov,alg,minFreq,alphabet ):
  """Return alternative allele only if different than ref and freq >= minFreq.
  """
  for base_count,base in zip( _get_base_counts( base_ref,alg,alphabet )[0],alphabet ):
    freq=base_count*1.0/cov
    if base!=base_ref and freq >= minFreq: return (base,freq) # base!=base_ref and 
  
def parse_mpileup( fpath,out_base,minDepth,minFreq,indels,alphabet='ACGT' ):
  """
  """
  #select input
  if not fpath: #stdin if not infile
    handle=sys.stdin
  else:         #file
    handle=open(fpath)

  #select output
  if out_base:    #write to file, if specified
    out1 = open( '%s.SNPs.%s.%s.txt' % (out_base,minDepth,minFreq),'w')
  else:           #write to stdout
    out1 = sys.stdout

  #write header
  out1.write( '#contig\tposition\tcoverage\tref\talt\tfreq\n')

  #process mpileup
  contigs=[]
  totCov={}; totLen={}; pContig=pPos=0
  for line in handle: 
    #line=line.strip()
    contig,pos,base_ref,cov,alg,quals=line[:-1].split('\t')#; contig,pos,base,cov,alg,quals
    pos,cov=int(pos),int(cov)
    #check whether same or next contig
    if pContig!=contig:
      pContig=contig
      contigs.append(contig)
      totCov[contig]=totLen[contig]=0
      pPos=0
    if indels and pos-pPos!=1:
      lineOut='%s\t%s-%s\tDEL\t\t\t\n' % ( contig,pPos,pos, )
      out1.write( lineOut )#; print lineOut,
    pPos=pos
    #add cov of the base
    totCov[contig]+=cov
    totLen[contig]+=1 #pPos=pos
    #check for heterozygous only if coverage > minDepth
    if cov<minDepth: continue
    alt_allele=get_alt_allele( base_ref,cov,alg,minFreq,alphabet )
    if alt_allele:
      base,freq = alt_allele
      lineOut='%s\t%s\t%s\t%s\t%s\t%1.4f\n' % ( contig,pos,cov,base_ref,base,freq )
      out1.write( lineOut )#; print lineOut,
  
  if out_base:    #close outfile, if opened
    out1.close() 

def main():

  usage = "usage: samtools mpileup -f ref.fa bam | %prog [options]" 
  parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True

  parser.add_option("-i", dest="fpath",
                    help="input file [stdin]")
  parser.add_option("-o", dest="outBase",
                    help="output base fname [%default]")
  parser.add_option("-d", dest="minDepth", default=10,  type=int,
                    help="minimal depth [%default]")
  parser.add_option("-f", dest="minFreq",  default=0.8, type=float,
                    help="min frequency of alternative base [%default]")
  parser.add_option("-n", dest="indels",   default=False, action="store_true", 
                    help="report indels [%default]")
  parser.add_option("-v", dest="verbose",  default=False, action="store_true" )
  
  ( o, fPaths ) = parser.parse_args()
  
  #parse pileup
  parse_mpileup( o.fpath,o.outBase,o.minDepth,o.minFreq,o.indels )
  
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
