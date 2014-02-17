#!/usr/bin/env python
###
#   Filter fasta file according to given filter.
# USAGE: 
# filter_fasta.py -iPL429-contigs.fa -l1000 > PL429-contigs.1000.fa
###

import getopt,numpy,os,sys
from commands import getoutput
from datetime import datetime
from Bio import SeqIO
from optparse import OptionParser

def filter_fasta( fname,minLength=0,upper=0,simpleHeader=False,window=50 ):
  """
  """
  for r in SeqIO.parse( open(fname),'fasta' ):
    if minLength and len(str(r.seq))<minLength: continue
    #
    j=0
    seq='%s\n' % r.seq[j:j+window]
    while j<=len(r.seq):
      j+=window
      seq+='%s\n' % r.seq[j:j+window]
    
    if upper: seq=seq.upper()
    if simpleHeader:  print '>%s\n%s' % ( r.id,         seq )
    else:             print '>%s\n%s' % ( r.description,seq )

def main( argv ):
  """
  """
  parser = OptionParser() #allow_interspersed_args=True
  
  parser.add_option("-i", "--fpath", dest="fpath", default=False,
                    help="mpileup file. if not specified, read from stdin (piping)")
  parser.add_option("-l", "--minLength", dest="minLength", default=0, type=int,
                    help="define minimal length of transcript", metavar="INT")
  parser.add_option("-s", "--simpleHeader", action="store_true", dest="simpleHeader", default=False,
                    help="simplify header f.e. '101751 length 5971 cvg_62.5_tip_0' > '101751' ")
  parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                    help="don't print status messages to stdout")
  parser.add_option("-u", "--upper", action="store_true", dest="upper", default=False, 
                    help="save sequence in upper case" )
  parser.add_option("-w", "--lineLenght", dest="lineLenght", default=70, type=int,
                    help="define max seq line length", metavar="INT")              
  parser.add_option("-r", "--replace", action="store_true", dest="replace", default=False,
                    help="overwrite output files", metavar="TXT")
  
  ( options, fPaths ) = parser.parse_args()
  if options.verbose: print options, fPaths   #fnames -> files of paired-ends reads in GERALD format (may be gzipped) ie. CDC37_6h_read1_qseq.txt.gz CDC37_6h_read2_qseq.txt.gz
  if not options.fpath and not options.replace: 
    sys.exit( " Parameter error: -i/--fpath has to be specified!" )
  
  filter_fasta( options.fpath,options.minLength,options.upper,options.simpleHeader,options.lineLenght )

if __name__=='__main__': 
  t0=datetime.now()
  main( sys.argv[1:] )
  dt=datetime.now()-t0
  #sys.stderr.write( "#Time elapsed: %s\n" % dt )
