#!/usr/bin/env python2
"""
Print stats about reads aligned to particular species.
Contigs/chromosomes have to be named as: spcode_contigcode


"""

import os, sys
from optparse import OptionParser
from datetime import datetime
from numpy    import mean, std

def idxstats2species( fpath,outfn ):
    """Prints number of reads aligned to species
    based on samtools idxstats.
    supercont1.26 1550015648 0
    supercont1.27 22892 0
    supercont1.28 2022244 0
    supercont1.29 347904 0
    * 0 0 0
    """
    #select input
    if not fpath: #stdin if not infile
        handle = sys.stdin
    else:         #file
        handle = open(fpath)

    #select output
    if outfn:    #write to file, if specified
        out1 = open( outfn,'w')
    else:           #write to stdout
        out1 = sys.stdout

    # process idxstats
    algTotal = 0
    sp2count = {}
    contigs, lengths, algs = [], [], []
    for l in handle:
        contig,length,alg = l.split()[:3]
        sp         = contig.split('_')[0]
        length,alg = int(length),int(alg)
        algTotal  += alg

        if contig=='*':
            continue

        if not sp in sp2count:
            sp2count[sp] = { "algs":    [],
                             "contigs": [],
                             "lengths": [] }
        sp2count[sp]["algs"].append( alg )
        sp2count[sp]["contigs"].append( contig )
        sp2count[sp]["lengths"].append( length )

    # process entries for all species
    for sp in sorted( sp2count.keys() ):
        contigs = sp2count[sp]["contigs"]
        # total assembly
        gSize = sum( sp2count[sp]["lengths"] )

        # aligned reads
        algs  = sum( sp2count[sp]["algs"] )

        # mean coverage
        meanCov = algs*1.0/gSize
        
        #print info
        lines = "%s\n%s contigs\nReads per base:\t%.2f\nAligned reads: %s [%.2f%s]\n" % ( sp,len(contigs),meanCov,algs,algs*100.0/algTotal,'%' )
        out1.write( lines )

    if outfn:
        out1.close()

def main():
    version = "%prog 1.0"
    usage = "usage: samtools idxstats bam | %prog [options]" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True

    parser.add_option("-i", dest="fpath",
                      help="input file [stdin]")
    parser.add_option("-o", dest="outfn",
                      help="output fname [stdout]")
    parser.add_option("-v", dest="verbose",  default=False, action="store_true" )
  
    ( o, fPaths ) = parser.parse_args()
  
    #parse pileup
    idxstats2species( o.fpath,o.outfn )
  
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
