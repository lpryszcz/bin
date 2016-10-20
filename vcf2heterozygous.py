#!/usr/bin/env python
desc="""Report only heterozygous mutations from vcf. 

"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 9/10/2012
"""

import argparse, os, pysam, sys
from datetime import datetime

def info2dict( info ):
    """Return dictionary of VCF format column.
    DP=56;VDB=0.0001;AF1=0.4999;CI95=0.5,0.5;DP4=9,11,5,13;MQ=11;FQ=8.65;PV4=0.33,0.089,1,0.0019
    """
    infodict = {}
    for kv in info.split(";"):
        if not kv:
            continue
        k,v=kv.split("=")
        if ',' in v:
            v = [ float(x) for x in v.split(",") ]
        else:
            v = float(v)
        infodict[k] = v
    return infodict

def vcf2heterozygous( vcf,freqth,maxfreqth,qualth,indels,verbose ):
    """Report heterozygous mutations.
    """
    if freqth and not maxfreqth:
        maxfreqth = 1-freqth
    for l in vcf:
        l = l.strip()
        if not l or l.startswith("#"):
            print l
            continue
        contig,pos,id,ref,alt,qual,flt,info,fmt = l.split("\t")[:9]
        #filter indels
        if not indels and info.startswith("INDEL"):
            continue
        #skip by quality 
        if float(qual)<qualth:
            continue
        infodict = info2dict( info )
        #check freq
        freq = 0
        if "DP4" in infodict:
            freqs = infodict["DP4"]
            freq  = sum(freqs[:2]) / sum(freqs)
        if not freqth <= freq <= maxfreqth:
            continue
        print l
        

def main():
    usage  = "usage: %(prog)s [options]"
    parser = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )

    parser.add_argument("-v", dest="verbose",  default=False, action="store_true")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument("-i", dest="vcf",      default=sys.stdin, type=file, 
                        help="vcf file            [%(default)s]" )
    parser.add_argument("-f", dest="freq",      default=0.3, type=float, 
                        help="base frequency      [%(default)s to 1-%(default)s]" )    
    parser.add_argument("-m", dest="maxfreq",   default=0, type=float, 
                        help="max frequency       [No filter or 1-freq if freq specified]" )    
    parser.add_argument("-q", dest="qual",      default=20.0, type=float, 
                        help="min SNP quality     [%(default)s]" )
    parser.add_argument("-n", dest="indels",    default=False, action="store_true",
                        help="report indels       [%(default)s]" )
   
    o = parser.parse_args()
    
    if o.verbose: 
        sys.stderr.write( "Options: %s\n" % str(o) )

    vcf2heterozygous( o.vcf,o.freq,o.maxfreq,o.qual,o.indels,o.verbose )
        
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )