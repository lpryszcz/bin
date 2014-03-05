#!/usr/bin/env python
desc="""Convert sam to fastq/fastq.
NOTE: sam has to be sorted by read name if paired output requested!
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Dublin, 28/06/2012
"""

import argparse, os, sys
from Bio      import Seq
from datetime import datetime

def sam2fasta( handle,outlist,last,first,paired ):
    """Return tuple of fastq and bool of sam entry.
    Bool is true if sam entry is correctly aligned."""
    pname = 0
    fqs   = []    
    for sam in handle:
        #skip empty lines and @ lines
        sam = sam.strip()
        if not sam or sam.startswith("@"):
            continue
        #for reference look at http://samtools.sourceforge.net/samtools.shtml
        qname,flag,tname,pos,mapq,cigar,mrnm,mpos,tlen,seq,qual = sam.split("\t")[:11]
        #get binary flag
        flag = int(flag)
        #if mapped to reverse strand 
        if flag & 16:
            #get reverse complement of seq
            seq   = str( Seq.Seq( seq ).reverse_complement() )
        #trim
        if last:
            seq = seq[:last]
        if first:
            seq = seq[first:]
        #get fasta
        fa = ">%s\n%s\n" % ( qname,seq )
        #report fasta
        if paired:
            if pname==qname:
                fqs.append( (fa,flag) )
                report( fqs,outlist )
                fqs = []
            else:
                fqs = [ (fa,flag) ]
            pname=qname
        else:
            outlist[0].write( fa ) 

def report( fqs,outlist ):
    """
    """
    for fq,flag in fqs:
        #paired and second in pair
        if flag & 128: 
            outlist[1].write( fq )
        else:
            outlist[0].write( fq )    
            
def sam2fastq( handle,outlist,last,first,paired ):
    """Return tuple of fastq and bool of sam entry.
    Bool is true if sam entry is correctly aligned."""
    pname = 0
    fqs   = []
    for sam in handle:
        #skip empty lines and @ lines
        sam = sam.strip()
        if not sam or sam.startswith("@"):
            continue
        #for reference look at http://samtools.sourceforge.net/samtools.shtml
        qname,flag,tname,pos,mapq,cigar,mrnm,mpos,tlen,seq,qual = sam.split("\t")[:11]
        #get binary flag
        flag = int(flag)
        #if mapped to reverse strand 
        if flag & 16:
            #get reverse complement of seq
            seq   = str( Seq.Seq( seq ).reverse_complement() )
            qual  = qual[::-1]
        #trim
        if last:
            seq,qual = seq[:last],qual[:last]
        if first:
            seq,qual = seq[first:],qual[first:]
        #get fastq
        fq = "@%s\n%s\n+\n%s\n" % ( qname,seq,qual )
        #report fastq
        if paired:
            if pname==qname:
                fqs.append( (fq,flag) )
                report( fqs,outlist )
                fqs = []
            else:
                fqs = [ (fq,flag) ]
            pname=qname
        else:
            outlist[0].write( fq )

def main():

    usage  = "samtools view BAM [region(s)] | %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-f", dest="fasta",   default=False, action="store_true",
                        help="output fasta       [fastq]")
    parser.add_argument("-o", dest="outbase", default=0, 
                        help="define output name [stdout]")
    parser.add_argument("-p", dest="paired",  default=False, action="store_true",
                        help="separate pairs     [%(default)s]")
    parser.add_argument("-l", dest="last",    default=0, type=int,
                        help="last base to keep 0-based  [entire read]")
    parser.add_argument("-s", dest="first",   default=0, type=int,
                        help="first base to keep 0-based [entire read]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #ext
    ext = "fastq"
    if o.fasta:
        ext = "fa"
        
    #define output files
    if o.paired:
        if not o.outbase:
            parser.error( "Specify output base name if paired output required!" )
        #define outnames
        fn1 = "%s_1.%s" % (o.outbase,ext)
        fn2 = "%s_2.%s" % (o.outbase,ext)
        #get opened files
        outlist = [ open(fn1,"w"),open(fn2,"w") ]
    elif o.outbase:
        fn1 = "%s.%s" % (o.outbase,ext)
        outlist = [ open(fn1,"w"),]
    else:
        outlist = [ sys.stdout,]
    
        
    if o.fasta:
        sam2fasta( sys.stdin,outlist,o.last,o.first,o.paired )
    else:
        sam2fastq( sys.stdin,outlist,o.last,o.first,o.paired )
        
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )    