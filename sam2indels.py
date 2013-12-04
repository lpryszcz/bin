#!/usr/bin/env python
desc="""Report SAM entries with  indel of size >= indelTh. 
NOTE: SAM has to be readname sorted. 
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 16/1/2013
"""

import argparse, os, re, sys
from Bio      import Seq
from datetime import datetime

#patM=re.compile(r'\d+M')
patIndel=re.compile(r'\d+I|\d+D')

def sam2indels( handle,indelTh,verbose ):
    """Reports sam algs having indel of size >= indelTh."""
    i=k=ptname=0
    if verbose:
        sys.stderr.write("Processing alignments...\n")
    for sam in handle:
        #skip empty lines and @ lines
        if sam.startswith("@"):
            sys.stdout.write(sam)

        i += 1
        #for reference look at http://samtools.sourceforge.net/samtools.shtml
        if len(sam.split('\t'))<11:
            continue
            
        qname,flag,tname,pos,mapq,cigar,mrnm,mpos,tlen,seq,qual = sam.split("\t")[:11]
        #get number of substitutions
        mtchs  = [ int(m[:-1]) for m in patIndel.findall( cigar ) ]
        #get fraction of read aligned
        #algfr = sum(mtchs)*1.0/len(seq)
        if mtchs and max(mtchs) >= indelTh:
            #report fastq
            #print cigar,mtchs
            sys.stdout.write( sam )
            k += 1
    
    sys.stderr.write( "%s alignments processed. %s [%.3f%s] having at least one indel of %s bases.\n" % (i,k,k*100.0/i,'%',indelTh) )
    
def main():

    usage  = "usage: samtools view BAM [region(s)] | %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="inputs",   nargs='+',default=[sys.stdin,], type=file, 
                        help="input file(s)           [stdin]")
    parser.add_argument("-s", dest="size",   default=10, type=int,
                        help="min. indel size in read [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    for handle in o.inputs:
        #if hasattr(handle,'name'):
        sys.stderr.write( "%s\t" % handle.name )
        sam2indels( handle,o.size,o.verbose )
        
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
    