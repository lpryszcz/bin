#!/usr/bin/env python
desc="""Report SAM entries with mismatches >= m.

"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 26/10/2012
"""

import argparse, os, re, sys
from Bio      import Seq
from datetime import datetime

patM=re.compile(r'\d+M')

def sam2mismatches_single( handle,out,algTh,verbose ):
    """Return tuple of fastq and bool of sam entry.
    Bool is true if sam entry is correctly aligned."""
    out = open( out+".fastq","w")
    i=k=ptname=0
    if verbose:
        sys.stderr.write("Processing alignments...\n")
    for sam in handle:
        #skip empty lines and @ lines
        if sam.startswith("@"):
            continue
        i += 1
        #for reference look at http://samtools.sourceforge.net/samtools.shtml
        qname,flag,tname,pos,mapq,cigar,mrnm,mpos,tlen,seq,qual = sam.split("\t")[:11]
        if verbose and ptname!=tname:
            ptname = tname
            sys.stderr.write( "  %s          \r" % tname )
        #get number of substitutions
        mtchs  = [ int(m[:-1]) for m in patM.findall( cigar ) ]
        #get fraction of read aligned
        algfr = sum(mtchs)*1.0/len(seq)
        if algfr < algTh:
            #if mapped to reverse strand 
            if int(flag) & 16:
                #get reverse complement of seq
                seq   = str( Seq.Seq( seq ).reverse_complement() )
                qual  = qual[::-1]
            #store fastq
            out.write( "@%s/1\n%s\n+\n%s\n" % ( qname,seq,qual ) )
            k += 1
    
    sys.stderr.write( "%s alignments processed.\n%s [%.3f%s] having <%.2f length aligned.\n" % (i,k,k*100.0/i,'%',algTh) )

def save_fastq( qname,out1,out2,flag,seq,qual,overlap ):
    """
    """
    #if mapped to reverse strand 
    if int(flag) & 16:
        #get reverse complement of seq
        seq   = str( Seq.Seq( seq ).reverse_complement() )
        qual  = qual[::-1]
    #split seq
    m = len(seq)/2
    o = int( round(overlap*len(seq)) )
    #store fastq
    out1.write( "@%s/1\n%s\n+\n%s\n" % ( qname,seq[:m+o],qual[:m+o] ) )
    out2.write( "@%s/2\n%s\n+\n%s\n" % ( qname,seq[m-o:],qual[m-o:] ) )
    
def sam2mismatches( handle,out,algTh,overlap,verbose ):
    """Return tuple of fastq and bool of sam entry.
    Bool is true if sam entry is correctly aligned."""
    out1 = open( out+"_1.fastq","w")
    out2 = open( out+"_2.fastq","w")
    i=k=ptname=0
    if verbose:
        sys.stderr.write("Processing alignments...\n")
    for sam in handle:
        #skip empty lines and @ lines
        if sam.startswith("@"):
            continue
        i += 1
        #for reference look at http://samtools.sourceforge.net/samtools.shtml
        qname,flag,tname,pos,mapq,cigar,mrnm,mpos,tlen,seq,qual = sam.split("\t")[:11]
        if verbose and ptname!=tname:
            ptname = tname
            sys.stderr.write( "  %s          \r" % tname )
        #get number of substitutions
        mtchs  = [ int(m[:-1]) for m in patM.findall( cigar ) ]
        #get fraction of read aligned
        algfr = sum(mtchs)*1.0/len(seq)
        if algfr < algTh:
            save_fastq( qname,out1,out2,flag,seq,qual,overlap )
            k += 1
    
    sys.stderr.write( "%s alignments processed.\n%s [%.3f%s] having <%.2f length aligned.\n" % (i,k,k*100.0/i,'%',algTh) )
    
def main():

    usage  = "usage: samtools view BAM [region(s)] | %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-o", dest="out",     default="",
                        help="define output basename         [mandatory]")
    parser.add_argument("-a", dest="algTh",   default=0.8, type=float,
                        help="max. fraction of read aligned  [%(default)s]")
    parser.add_argument("-b", dest="overlap", default=0.1, type=float,
                        help="franction of split read overlapping [%(default)s]")
    parser.add_argument("-s", dest="single",  default=False, action="store_true",
                        help="report single reads            [%(default)s]")    
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    if not o.out:
        parser.error( "Output name has to be specified!" )

    if o.single:
        sam2mismatches_single( sys.stdin,o.out,o.algTh,o.verbose )
    else:
        sam2mismatches( sys.stdin,o.out,o.algTh,o.overlap,o.verbose )
        
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
    