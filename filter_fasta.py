#!/usr/bin/env python
desc="""Filter fasta file according to given filter.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 2012
"""

import argparse, numpy, os, sys
from datetime import datetime
from Bio import SeqIO

def filter_fasta(handle, minLength=0, upper=0, simpleHeader=False, window=50, \
                 skip_xs=0, skip_ns=0, totalSize=0):
  """
  """
  totsize = 0
  for r in SeqIO.parse(handle, 'fasta'):
    if minLength and len(r)<minLength: continue
    #
    j=0
    seq='%s\n' % r.seq[j:j+window]
    while j<=len(r.seq):
      j+=window
      seq+='%s\n' % r.seq[j:j+window]
    #skip if too many Xs / Ns in sequence
    if skip_xs and 1.0 * r.seq.count('X') / len(r) > skip_xs:
        continue
    if skip_ns and 1.0 * r.seq.count('N') / len(r) > skip_ns:
        continue
    #
    if upper: seq=seq.upper()
    if simpleHeader:  print '>%s\n%s' % ( r.id,         seq )
    else:             print '>%s\n%s' % ( r.description,seq )
    # break if total size exceeded
    totsize += len(r)
    if totalSize and totsize>totalSize:
      break

def main():

    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-i", "--input", default=sys.stdin,type=file, 
                        help="input stream [stdin]")
    parser.add_argument("-l", "--minLength", dest="minLength", default=0, type=int,
                        help="define minimal length of transcript")
    parser.add_argument("-t", "--totalSize", default=0, type=int,
                        help="max size to output")
    parser.add_argument("-s", "--simpleHeader", action="store_true", default=False,
                        help="simplify header f.e. '101751 length 5971 cvg_62.5_tip_0' > '101751' ")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="don't print status messages to stdout")
    parser.add_argument("-u", "--upper", action="store_true", default=False, 
                        help="save sequence in upper case" )
    parser.add_argument("-w", "--lineLenght", default=70, type=int,
                        help="define max seq line length")              
    parser.add_argument("-r", "--replace", action="store_true",  default=False,
                        help="overwrite output files")
    parser.add_argument("-x", "--skip_xs", type=float, default=0,
                        help="skip sequences with Xs")
    parser.add_argument("-n", "--skip_ns", type=float, default=0,
                        help="skip sequences with ns")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
  
    filter_fasta(o.input, o.minLength, o.upper, o.simpleHeader, o.lineLenght, \
                 o.skip_xs, o.skip_ns, o.totalSize)

if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("\nI/O error({0}): {1}\n".format(e.errno, e.strerror))  
    dt  = datetime.now()-t0
    #sys.stderr.write("#Time elapsed: %s\n" % dt)
    
