#!/usr/bin/env python2
"""Compute best reciprocal blast hits for two multifasta protein files.
Build blast db and launch blastp in 2 directions. 
Keep only best reciprocal matches.

USAGE:
python BRH.py f1.faa f2.faa

DEPENDENCIES:
blastall
Biopython

Author:
lpryszcz@crg.es

!!! Deprecated, use brh.py instead !!!
"""

from Bio import SeqIO
import os, sys
from datetime import datetime

def parse_blastout( fn ):
  q2t = {}
  for line in open( fn ):
    #CPAG_00001.1	orf19.3540	49.80	249	87	3	46	262	13	255	3e-39	 156
    q,t,identity,alength,mismatches,gaps,qstart,qend,sstart,send,evalue,score=line[:-1].split('\t')
    if q in q2t:
      continue
    q2t[q]=(t,identity,alength,mismatches,gaps,qstart,qend,sstart,send,evalue,score)
  return q2t

def main():
  if len( sys.argv )< 3:
    sys.exit( "Provide two multifasta files as input." )
  
  fn1,fn2 = sys.argv[1:3]
  
  for fn in (fn1,fn2):
    if not os.path.isfile( fn ):
      sys.exit( "No such file: %s" % fn )

  #build db
  sys.stderr.write( "Building databases...\n" )
  for fn in (fn1,fn2):
    if not os.path.isfile( "%s.pin" % fn ): #CANAL.proteins.fa.pin
      os.system( "formatdb -oT -pT -i %s" % fn )
      
  #run blastp
  sys.stderr.write( "Executing blastp...\n" )
  for f1,f2 in zip( (fn1,fn2),(fn2,fn1) ):
    if not os.path.isfile( "%s.%s.blastout.txt" % ( f1,f2 ) ):
      os.system( "blastall -pblastp -e1e-05 -m8 -a3 -i %s -d %s > %s.%s.blastout.txt" % ( f1,f2,f1,f2 ) )
      
  #compare matches
  q2t1 = parse_blastout( "%s.%s.blastout.txt" % ( fn1,fn2 ) )
  q2t2 = parse_blastout( "%s.%s.blastout.txt" % ( fn2,fn1 ) )
  
  for q in q2t1:
    t=q2t1[q][0] #target,identity,alength,mismatches,gaps,qstart,qend,sstart,send,evalue,score
    #if t in other list as q
    if t in q2t2:
      #if BRH
      if q2t2[t][0]==q:
        sys.stdout.write( "%s\t%s\t%s\t%s\t%s\t%s\n" % (q,t,q2t1[q][1],'\t'.join(q2t1[q][-2:]),q2t2[t][1],'\t'.join(q2t2[t][-2:]) ) )

def brh( fn1, fn2 ):
  """
  """
  for fn in (fn1,fn2):
    if not os.path.isfile( fn ):
      sys.exit( "No such file: %s" % fn )

  #build db
  sys.stderr.write( "Building databases...\n" )
  for fn in (fn1,fn2):
    if not os.path.isfile( "%s.pin" % fn ): #CANAL.proteins.fa.pin
      os.system( "formatdb -oT -pT -i %s" % fn )
      
  #run blastp
  sys.stderr.write( "Executing blastp...\n" )
  for f1,f2 in zip( (fn1,fn2),(fn2,fn1) ):
    if not os.path.isfile( "%s.%s.blastout.txt" % ( f1,f2 ) ):
      os.system( "blastall -pblastp -e1e-05 -m8 -a3 -i %s -d %s > %s.%s.blastout.txt" % ( f1,f2,f1,f2 ) )
      
  #compare matches
  q2t1 = parse_blastout( "%s.%s.blastout.txt" % ( fn1,fn2 ) )
  q2t2 = parse_blastout( "%s.%s.blastout.txt" % ( fn2,fn1 ) )
  
  for q in q2t1:
    t=q2t1[q][0] #target,identity,alength,mismatches,gaps,qstart,qend,sstart,send,evalue,score
    #if t in other list as q
    if t in q2t2:
      #if BRH
      if q2t2[t][0]==q:
        sys.stdout.write( "%s\t%s\t%s\t%s\t%s\t%s\n" % (q,t,q2t1[q][1],'\t'.join(q2t1[q][-2:]),q2t2[t][1],'\t'.join(q2t2[t][-2:]) ) )
        
def main1():
  usage   = "usage: %prog [options] [ > brh.txt ]" #arg1 arg2
  version = "%prog 1.0"
  parser  = OptionParser( usage=usage,version=version ) #allow_interspersed_args=True
  
  parser.add_option( "-i", dest="infile",         default="",  type=str, 
                     help="file1 with ordered genes     [mandatory]")
  parser.add_option( "-j", dest="infile2",        default="",  type=str, 
                     help="file2 with ordered genes. If not specified, within species synteny is calculated")
  #parser.add_option( "-d", dest="dbdir",          default="blastdb",  type=str, 
  #                   help="directory to store blast db    [%default ]")

  parser.add_option( "-m", dest="minCommon",      default=2,   type=int, 
                     help="min common loci                [%default    ]")
  parser.add_option( "-s", dest="sharedFaction",  default=0.01, type=float, 
                     help="min fraction of shared loci    [%default ]" )
  parser.add_option( "-t", dest="paired",         action="store_true", default=False,
                     help="ignore singleton genes         [%default]")
  parser.add_option( "-v", dest="verbose",        action="store_true", default=False )
                          
  algOpt = OptionGroup( parser, "Alignment options [not implemented]", )
  algOpt.add_option( "-n", dest="match",          default= 5, type=int, help="[ %default]" )
  algOpt.add_option( "-o", dest="mismatch",       default=-1, type=int, help="[%default]" )
  algOpt.add_option( "-p", dest="gapOpen",        default=-5, type=int, help="[%default]" )
  algOpt.add_option( "-r", dest="gapExtend",      default=-2, type=int, help="[%default]" )
  parser.add_option_group( algOpt )
  
  ( o, args ) = parser.parse_args()
  
  #check if file(s) exist
  if not o.infile:
    parser.error( "Specify input file" )
  if not os.path.isfile( o.infile ): 
    parser.error( "No such file: %s" % o.infile )
  if o.infile2 and not os.path.isfile( o.infile2 ): 
    parser.error( "No such file: %s" % o.infile2 )
  if o.conversion and not os.path.isfile( o.conversion ): 
    parser.error( "No such file: %s" % o.conversion ) 

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
