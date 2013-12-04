#!/usr/bin/env python
"""Split bam file into two, for + and - strand transcripts.
BAM has to be sorted.

USAGE:
bam2stranded.py bam

Transcript from 
+ strand are:
 83 163
- strand:
 99 141
"""

import os,sys
import subprocess

def bam2stranded( bam ):
  """ """
  #get stdin
  print "Reading from:"
  args = [ 'samtools','view','-h','-F 4',bam ]; print ""," ".join( args )
  proc = subprocess.Popen( args,stdout=subprocess.PIPE )
  #define stdout
  print "Writing to:"
  args1 = [ 'samtools','view','-Sb','-o%s.+.bam' % bam,'-' ]; print " +"," ".join( args1 )
  proc1 = subprocess.Popen( args1,stdin=subprocess.PIPE )
  args2 = [ 'samtools','view','-Sb','-o%s.-.bam' % bam,'-' ]; print " -"," ".join( args2 )
  proc2 = subprocess.Popen( args2,stdin=subprocess.PIPE )
  out   = [ proc1.stdin,proc2.stdin ]
  #process stdin
  print "Processing chromosomes..."
  pContig=""
  for line in proc.stdout:
    #write header lines to both outputs
    if line.startswith('@'):
      out[0].write( line )
      out[1].write( line )
      continue
    #read flag and chromosome
    read,flag,contig = line.split("\t")[:3]
    if contig!=pContig:
      pContig=contig
      sys.stderr.write( " %s     \r" % contig )
    #get binary flag
    bflag = bin(int(flag))
    #find from which strand transcript comes from
    if len(bflag)>8 and bflag[-7]=='1': #first in pair
      if bflag[-5]=='1': #mapped to -
        out[0].write( line )
      else: #mapped to +
        out[1].write( line )
    #second in pair
    elif len(bflag)>6 and bflag[-5]=='1': #mapped to -
      out[1].write( line )
    else: #mapped to +
      out[0].write( line )
  
  #close subprocesses
  out[0].close()
  out[1].close()
  #proc1.terminate()
  #proc2.terminate()
  
  print "Indexing..."
  os.system( 'samtools index %s.+.bam' % bam )
  os.system( 'samtools index %s.-.bam' % bam )
  print "Done!"
  
if __name__=="__main__":
  cmd = "USAGE: bam2stranded.py bam"
  fnames = sys.argv[1:]
  if not fnames:
    sys.exit( cmd )
  
  for bam in fnames:
    if not os.path.isfile( bam ):
      sys.exit( cmd + "\nNo such file: %s" % bam )
      
    bam2stranded( bam )
  
  

      

