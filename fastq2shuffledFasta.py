#!/usr/bin/env python2
"""
for f in *_R1.fastq.gz; do s=`echo $f | cut -d"_" -f1`; echo `date` $s; fastq2shuffledFasta.py -l31 -q20 > ../assemblies/identity_confirmation/mapsembler/$s.fa ${s}_*.gz 2> filterReads.log; done
"""
 
import gzip, os, sys
from Bio      import Seq
from datetime import datetime
from optparse import OptionParser
  
def qseq2fasta( file,minLen,qualityTh,reverse,ASCII_offset=33,sep='N' ):
  """Process each line of GERALD output and Return FastQ str. Cut seq @ first '.' position.
  Also, check for quality if qualityTh defined.
  Return None if length of the seq < minLen or an error occured during line processing.
  """
  name  = file.next()[:-1]
  seq   = file.next()[:-1]
  sep   = file.next()[:-1]
  quals = file.next()[:-1]
  
  #format name - @HWI-ST227:145:C06RAACXX:7:1101:1156:2148 2:Y:0:AACT > @HWI-ST227:145:C06RAACXX:7:1101:1156:2148/2
  #if len( name.split() ) > 1:
  #  name = '%s/%s' % ( name.split()[0],name.split()[1][0] )
  
  #clip seq & quals @ . ( unknown base )
  if sep in seq:
    pos=seq.index(sep)
    seq,quals=seq[:pos],quals[:pos]
    if len(seq)<minLen or not seq: 
      return
  
  #convert PHRED64 to PHRED33
  #quals=''.join( [ chr(ord(q)-31) for q in quals ] )

  #cut sequence & quals @ quality
  if qualityTh:
    pos=0
    for q in quals:
      phredQ=ord(q)-33 #PHRED+33 encoding
      if phredQ<qualityTh: 
        seq,quals=seq[:pos],quals[:pos]
        if len(seq)<minLen or not seq:  
          return
        break
      pos+=1
  
  if reverse:
    seq = str(Seq.Seq(seq).reverse_complement())
  #return fasta
  return seq #'>%s\n%s\n' % ( name,seq ) #'%s\n%s\n+\n%s\n' % ( name,seq,quals )
  
def qseq2shuffledFasta( fPathF,fPathR,minLen,qualityTh,reverse,ASCII_offset=33,gzlipEndings=('.gz',) ):
  """Convert GERALD files (fPathF and fPathR) to FastQ (outFileF,outFileR) with quality and read length filtering.
  Trim @ first base with quality < qualityTh, distard reads with lenth < minLen. 
  Write singletons FastQ to outUnpaired if exists.
  Write filtered pairs to combinedOutFile if exists.
  """
  #open input files
  i=j=pI=0
  if fPathF.endswith(gzlipEndings) and fPathR.endswith(gzlipEndings):
    fileF=gzip.open( fPathF,'rb' )
    fileR=gzip.open( fPathR,'rb' )
  else:
    fileF=open( fPathF,'rb' )
    fileR=open( fPathR,'rb' )
    
  #parse files
  step=10**4
  try:
    while 1:
      i+=1
      rec1 = qseq2fasta( fileF,minLen,qualityTh,reverse,ASCII_offset )
      rec2 = qseq2fasta( fileR,minLen,qualityTh,reverse,ASCII_offset )
      if rec1 and rec2:
        j+=1
        sys.stdout.write( '>%s/1\n%s\n>%s/2\n%s\n' % (i,rec1,i,rec2) )
        
      #print info
      if i>pI:
        pI+=step
        sys.stderr.write( "%s done. %s successes [%.2f%s]\r" % (i,j,j*100.0/i,'%') )
        
  except StopIteration, e: 
    sys.stderr.write( "Processed pairs: %s. Successfull: %s [%.2f%s].\n" % ( i,j,j*100.0/i,'%' ) )

def main():
  usage  = "usage: %prog [options] fq1 fq2"
  parser = OptionParser(usage=usage) #allow_interspersed_args=True
  
  parser.add_option("-l", dest="minLen", default=31, type=int,
                    help="min read lenght (shorter after quality trimming are removed) [default: %default]" )
  parser.add_option("-q", dest="qualityTh", default=0, type=int,
                    help="read is clipped @ first base having PHRED quality lower than [default: %default]" )
  parser.add_option("-r", dest="reverse",   default=False, action='store_true', 
                    help="reverse complement [default: %default]" )
  #parser.add_option("-t", dest="PHRED", default=64,
  #                  help="select input ASCII offset: Sanger -> PHRED33, Solexa -> PHRED64 [default: %default]")     
  
  ( o, fPaths ) = parser.parse_args()
  sys.stderr.write( "Options: %s\nFiles to be processed: %s\n" % ( o, ', '.join( fPaths ) ) )
  
  ###check input parameters
  #file input are files
  for fpath in fPaths:
    if not os.path.isfile( fpath ): 
      sys.exit( "Error! No such file: %s" % fpath )
  
  #process pairs of files
  for fPathF,fPathR in zip( fPaths[::2],fPaths[1::2] ):
    qseq2shuffledFasta( fPathF,fPathR,o.minLen,o.qualityTh,o.reverse )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
