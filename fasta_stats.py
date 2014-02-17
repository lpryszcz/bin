#!/usr/bin/env python
###
# Prints info about multifasta file (dna):
# 1.number of contigs
# 2.size
# 3.GC
# 4.Ns
# USAGE: fasta_stats.py *.scafSeq.1000.fa
#
import os, sys
from Bio import SeqIO
from datetime import datetime

def fasta_stats( fn ):
  statsFn=fn+'.stats'
  #return content of .stats file if exists and younger than fasta
  if os.path.isfile(statsFn) and os.stat(fn).st_mtime < os.stat(statsFn).st_mtime:
    line=open(statsFn).readline()
    line=line.strip()
    line='%s\t' % fn + '\t'.join( line.split('\t')[1:] )
    return line
  #if not, generate that file
  lengths=[]; lengths1000=[]
  contigs=contigs1000=baseErrors=0
  #count bases frequencies
  bases={'A':0,'C':0,'G':0,'T':0,'N':0}
  for r in SeqIO.parse( open(fn),'fasta' ):
    contigs+=1
    seq=str(r.seq)
    seq=seq.upper()
    lengths.append( len(seq) )
    if len(seq)>1000: 
      contigs1000+=1
      lengths1000.append( len(seq) )
    for base in seq:
      try:    bases[base]+=1
      except: baseErrors+=1
  if not lengths:
    return fn+'\tError: No sequences!'
  #calculate GC
  if bases['A']+bases['T']: GC=(bases['G']+bases['C'])*100.0/(bases['A']+bases['C']+bases['G']+bases['T'])
  else:                     GC=0
  
  #N50 & N90
  size=sum(lengths) #sum(bases.itervalues())
  lengthSum=0
  n50=n90=0
  lengths.sort( reverse=True )
  for l in lengths:
    lengthSum += l
    if not n50 and lengthSum>=0.5*size:
      n50=l
    if lengthSum>=0.9*size:
      n90=l 
      break
      
  #print output
  line='%s\t%s\t%s\t%.3f\t%s\t%s\t%s\t%s\t%s\t%s' % ( fn,contigs,size,GC,contigs1000,sum(lengths1000),n50,n90,bases['N'],lengths[0] )
  try:
  	out=open(statsFn,'wb'); out.write( line+'\n' ); out.close()
  except IOError as e:
  	sys.stderr.write("%s\n" % e)
  return line

def main( fnames ):
  if not fnames:
    sys.exit( "USAGE: \nfasta_stats.py *.contigs.fa" )
    
  print '#fn\tcontigs\tbases\tGC [%]\tcontigs >1kb\tbases in contigs >1kb\tN50\tN90\tNs\tlongest'
  for fn in fnames:
    print fasta_stats( fn )
  
if __name__=='__main__': 
  T0=datetime.now()
  main( sys.argv[1:] )
  print "#Elapsed time: %s" % ( datetime.now()-T0, )
