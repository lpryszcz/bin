#!/usr/bin/env python
desc="""Concatenate multiple fasta files into phylip alignment.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Girona/Cracow, 6/11/2012
"""

import argparse, os, sys
from datetime import datetime
from Bio      import SeqIO, SeqRecord

def fasta2concatenated(out, files, splitname, verbose):
    """Concatenate multiple fasta entries into one.
    Concatenate multiple fasta files into one.
    """
    if verbose:
        sys.stderr.write("Concatenating fastas...\n")

    records = []
    for f in files:
        if verbose:
            sys.stderr.write(" %s          \r" %f.name )
        #define seq name
        name = f.name
        if splitname:
            name = f.name.split(".")[0]
        #get all entries from fasta file
        seq = ""
        for r in SeqIO.parse(f, "fasta"):
            seq += r.seq
        #create sequence record
        sr = SeqRecord.SeqRecord(seq, id=name) 
        records.append(sr)
    #store all SeqRecords
    SeqIO.write(records, out, "phylip")

def main():

    usage  = "usage: %(prog)s [options] -i fasta1 [fasta2 ... fastaN]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="files",   nargs="+", type=file, required=True, 
                        help="fasta fnames        [%(default)s]") 
    parser.add_argument("-o", dest="out",     default=sys.stdout, type=argparse.FileType('w'), 
                        help="define output name  [stdout]")
    parser.add_argument("-s", dest="splitname", default=False, action="store_true",
                        help="split seq name at dot [%(default)s]")    
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #go
    fasta2concatenated( o.out,o.files,o.splitname,o.verbose )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
  
