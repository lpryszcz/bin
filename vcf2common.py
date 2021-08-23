#!/usr/bin/env python2
desc="""Get common SNPs
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Sosnowiec, 8/11/2012
"""

import argparse, os, sys
from datetime import datetime
from Bio      import SeqIO

def update_variants( s,vcf,variants,verbose ):
    """Return dictionary of variants
    with 0-based, half-open coordinates (BED,python).
    VCF is 1-based (GFF like) ?
    """
    for l in vcf:
        l = l.strip()
        if not l or l.startswith("#"):
            continue
        lData = l.split("\t")
        contig,pos,id,ref,alt,qual,flt = lData[:7]
        #define variant descriptor
        v="%s_%s_%s_%s"%(contig,pos,ref,alt)
        if len(lData)>9:
            v += "_"+"_".join(lData[10:])
        #add variants
        if v not in variants:
            variants[v] = [s]
        else:
            variants[v].append(s)
    return variants
    
def vcf2common( out,vcfs,minmt,splitname,unique,verbose ):
    """
    """
    #sample2variants = {}
    variants = {}
    if verbose:
        sys.stderr.write( "Loading variants...\n" )
    #load vcf file(s)
    for vcf in vcfs:
        s = vcf.name
        if splitname:
            s = s.split(".")[0]
        if verbose:
            sys.stderr.write( " %s -> %s        \n" % (vcf.name,s) )
        #load vcf
        variants = update_variants( s,vcf,variants,verbose )
        #sys.stderr.write('   %s\n'%len(variants))
    #report unique or 
    for v,samples in variants.items():
        if unique and len(samples)==1 or not unique and len(samples)>=minmt:
            line = "%s\t%s\n" % (" ".join(samples),"\t".join(v.split("_")))
            sys.stdout.write(line)

        
def main():

    usage  = "usage: %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="vcfs",    nargs="+", type=file, 
                        help="vcf file(s)         [mandatory]") 
    parser.add_argument("-o", dest="out",     default=sys.stdout, type=argparse.FileType('w'), 
                        help="define output name  [stdout]")
    parser.add_argument("-n", dest="minmt",   default=1, type=int, 
                        help="min. common SNPs    [%(default)s]")
    parser.add_argument("-s", dest="splitname", default=False, action="store_true",
                        help="split name at dot   [%(default)s]")    
    parser.add_argument("-u", dest="unique",  default=False, action="store_true",
                        help="report unique       [%(default)s]")    
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #go
    vcf2common( o.out,o.vcfs,o.minmt,o.splitname,o.unique,o.verbose )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
  



