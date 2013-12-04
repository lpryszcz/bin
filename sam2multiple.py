#!/usr/bin/env python
desc="""Report SAM entries with multiple matches. 
NOTE: SAM has to be readname sorted. 
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 16/1/2013
"""

import argparse, os, sys
from datetime import datetime

def sam2multiple( handle,verbose ):
    """Report multiple sam algs for given read."""
    i = k = 0
    pqname = psam = ""
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
            
        #report sam if same name as previous
        if qname == pqname and ptname == tname:
            if psam:
                k += 1
                sys.stdout.write( psam )
            sys.stdout.write( sam )
            psam = ""
        #store psam only if new qname
        else:
            psam = sam
            
        #store qname
        pqname = qname
        ptname = tname
    
    sys.stderr.write( "%s alignments processed. %s [%.3f%s] having multiple matches.\n" % (i,k,k*100.0/i,'%') )
    
def main():

    usage  = "usage: samtools view BAM [region(s)] | %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="inputs",   nargs='+',default=[sys.stdin,], type=file, 
                        help="input file(s)           [stdin]")
    #parser.add_argument("-n", dest="size",   default=2, type=int,
    #                    help="min. indel size in read [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    for handle in o.inputs:
        #if hasattr(handle,'name'):
        sys.stderr.write( "%s\t" % handle.name )
        sam2multiple( handle,o.verbose )
        
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
    