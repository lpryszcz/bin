#!/usr/bin/env python
"""
Look for telomeric sequence in genomic reads from Illumina.
USAGE:
out=telomers.txt; 
for fn in */q10_1.fastq; do 
  d=`echo $fn | cut -d'/' -f1`;
  echo $d >> $out; 
  echo "find_telomers.py $d/q10*merFreq.txt"; 
  find_telomers.py $d/q10*merFreq.txt >> $out; 
done
"""

import sys
import re
from copy     import copy
from Bio.Seq  import Seq
from datetime import datetime

def get_kmers( seq,k=5,basesInMerTh=2 ):
  """
  """
  mers=set()
  for i in range(len(seq)-k):
    mer=seq[i:i+k]
    #check if mer is homozygous, skip if so
    bases=set( base for base in mer )
    if seq.count(mer)>1 and len(bases)>basesInMerTh: 
      mers.add(mer)
  return mers
  
def get_repeats( mer,seqs,minLen=7 ):
  """
  """
  #define pattern
  pat=re.compile(mer)
  repeats=set()
  #iterate through sequences
  for sid in seqs:
    for seq in seqs[sid]:
      pSpan=[]
      #store info about repeat if at least 2 matches
      for m in pat.finditer(seq):
        if pSpan:
          start,end=pSpan[0],m.span()[0]
          if end-start>=minLen:
            repeat=seq[ start:end ]
            #add repeat info
            repeats.add(repeat)
        pSpan=m.span()
  return repeats

def _unique( repeat,repeats ):
  """Return False if any repositioned repeat or
  its reverse complement in repeats.
  Otherwise return True.
  """
  for i in range( len(repeat) ):
    repeat=repeat[1:]+repeat[0]
    rcRepeat=str( Seq(repeat).reverse_complement() )
    if repeat in repeats or rcRepeat in repeats: return False
  return True

def get_unique_repeats( repeats ):
  """Return set of non-redundant repeats.
  """
  nrRepeats=set()
  for repeat in repeats:
    if _unique( repeat,nrRepeats ): nrRepeats.add( repeat )
  return nrRepeats
  
def count_repeats( seqs, repeats ):
  """
  """
  #print len(repeats),repeats
  #get only non-reduntant repeats
  repeats=get_unique_repeats( repeats )
  #print len(repeats),repeats
  
  #calculate stats
  matches={}
  for repeat in repeats:
    matches[repeat]=[0,0,0]
    for sid in seqs:
      #read pattern occurences from sequence id
      occurences=int( sid.split('-')[-1] )
      for seq in seqs[sid]:
        for i in range( 0,len(repeat),1 ):
          #create repeat+flanking 5bp
          repeat3=3*repeat
          _repeat=repeat[i:]+repeat[:i]+repeat3[i:i+len(seq)-len(repeat)]
          #print _repeat
          if _repeat in seq: 
            matches[ repeat ][0] += occurences
            matches[ repeat ][1] += 1
            break
          
  #add scores
  for repeat in matches:
    o,i,s=matches[repeat]
    s=o*len(repeat)
    matches[repeat]=o,i,s
    
  return matches  

def main():  
  k=5
  seqs={}
  for fn in sys.argv[1:]:
    f=open(fn)
    while 1:
      #read sequence id (contain occurencies info), sequence and it's reverse complement
      sid  = f.readline()[1:-1] #>5-1870
      seq  = f.readline()[:-1]  #GTTCAGCAGGAATGCCGAGATCGGAAGAGCGGTTCAGCAGGAA
      rSeq = f.readline()[:-1]  #TTCCTGCTGAACCGCTCTTCCGATCTCGGCATTCCTGCTGAAC
      #stop reading if not sequence id
      if not sid: break
      #add fn to sequnce id and store
      sid = '%s-%s' % (fn,sid)
      seqs[sid] = ( seq,rSeq )
    f.close()

  seqsOrg     = copy(seqs)
  sids        = seqs.keys()
  merProcessed= set()
  repeats     = set()
  for sid in sids:
    seq,rSeq = seqs.pop(sid)
    #print sid,seq,rSeq
    mers=get_kmers( seq,k )
    for mer in mers:
      #check if kmer already processed
      if mer in merProcessed: continue
      merProcessed.add( mer )
      #get repeats flanked by given kmer and add them to repeats
      for repeat in get_repeats( mer,seqs ): repeats.add( repeat )

  #here it will be cool to collapse overlapping repeats
  matches=count_repeats( seqsOrg, repeats )

  #print best repeats
  i=0
  #print '#reads\tunique\tscore\trepeat'
  for repeat in sorted( matches.keys(), key=lambda x: matches[x][2], reverse=True ):
    o,i,s=matches[repeat]
    if not int(i): break
    print '%s\t%s\t%s\t%s' % ( o,i,s,repeat, )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "  #Time elapsed: %s\n" % dt )
