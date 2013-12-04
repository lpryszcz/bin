#!/usr/bin/env python
"""
"""
 
import gzip, os, sys
from datetime import datetime
from optparse import OptionParser
  
def getFastQ( file,minLen,qualityTh,ASCII_offset_33=True,pair='',sep='N' ):
  """Process each line of GERALD output and Return FastQ str. Cut seq @ first '.' position.
  Also, check for quality if qualityTh defined.
  Return None if length of the seq < minLen or an error occured during line processing.
  """
  name  = file.next()[:-1]
  seq   = file.next()[:-1]
  sepli = file.next()[:-1]
  quals = file.next()[:-1]
  
  #format name - @HWI-ST227:145:C06RAACXX:7:1101:1156:2148 2:Y:0:AACT > @HWI-ST227:145:C06RAACXX:7:1101:1156:2148/2
  if len( name.split() ) > 1:
    name = '%s%s' % ( name.split()[0],pair )
  
  #clip seq & quals @ . ( unknown base )
  if minLen and sep in seq:
    pos = seq.index(sep)
    if pos<minLen:
      return
    seq,quals=seq[:pos],quals[:pos]
  
  #convert PHRED64 to PHRED33
  if not ASCII_offset_33:
    quals=''.join( [ chr(ord(q)-31) for q in quals ] )

  #cut sequence & quals @ quality
  if qualityTh:
    pos=0
    for q in quals:
      phredQ=ord(q)-33 #PHRED+33 encoding
      if phredQ<qualityTh: 
        seq,quals=seq[:pos],quals[:pos]
        break
      pos+=1
     
  if len(seq)<minLen or not seq:  
    return
  
  #return fastq
  return '%s\n%s\n+\n%s\n' % ( name,seq,quals )
  
def qseq2shuffledFastQ( fPathF,fPathR,minLen,qualityTh,ASCII_offset_33=True,upto=0,verbose=False,gzlipEndings=('.gz',) ):
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
      if upto and i>upto:
        break
      i+=1
      rec1 = getFastQ( fileF,minLen,qualityTh,ASCII_offset_33,'/1' )
      rec2 = getFastQ( fileR,minLen,qualityTh,ASCII_offset_33,'/2' )
      if rec1 and rec2:
        j+=1
        sys.stdout.write( '%s%s' % (rec1,rec2) )
        
      #print info
      if i>pI:
        pI+=step
        if verbose:
          sys.stderr.write( "%s done. %s successes [%.2f%s]\r" % (i,j,j*100.0/i,'%') )
        
  except StopIteration, e: 
    if verbose:
      sys.stderr.write( "Processed pairs: %s. Successfull: %s [%.2f%s].\n" % ( i,j,j*100.0/i,'%' ) )

def main():
  #usage  =
  parser = OptionParser( ) #allow_interspersed_args=True usage=usage 
  
  parser.add_option("-v", action="store_true", dest="verbose", default=False,
                    help="print status messages to stdout [default: %default]")
  parser.add_option("-l", dest="minLen", default=31, type=int,
                    help="min read lenght (shorter after quality trimming are removed) [default: %default]" )
  parser.add_option("-q", dest="qualityTh", default=3, type=int,
                    help="read is clipped @ first base having PHRED quality lower than [default: %default]" )
  parser.add_option("-t", dest="ASCII_offset_33", action="store_false", default=True,
                    help="use illumina/solexa quality encoding (ASCII offset of 64) [default: Sanger]")
  parser.add_option("-u", dest="upto", default=0, type=int,
                    help="process up to pairs [all]")
  
  ( o, fPaths ) = parser.parse_args()
  sys.stderr.write( "Options: %s\nFiles to be processed: %s\n" % ( o, ', '.join( fPaths ) ) )
  
  ###check input parameters
  #file input are files
  for fpath in fPaths:
    if not os.path.isfile( fpath ): 
      sys.exit( "Error! No such file: %s" % fpath )
  
  #process pairs of files
  for fPathF,fPathR in zip( fPaths[::2],fPaths[1::2] ):
    qseq2shuffledFastQ( fPathF,fPathR,o.minLen,o.qualityTh,o.ASCII_offset_33,o.upto,o.verbose )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
