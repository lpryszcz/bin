#!/usr/bin/env python
"""
Calculate mean coverage for given genome.
"""

import os, sys
from optparse import OptionParser
from datetime import datetime
from numpy    import mean, std

def coverage( fpath,outfn ):
    """Prints out coverage info base on samtools idxstats.
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
    
    contigs, lengths, algs = [], [], []
    for l in handle:
        contig,length,alg = l.split()[:3]
        length,alg = int(length),int(alg)
        if contig=='*':
            continue
        contigs.append( contig )
        lengths.append( length )
        algs.append( alg )

    #print stats
    meanCov = sum( algs )*1.0/sum( lengths )
    #contigsCov
    ctgsCov = [ a*1.0/l for l,a in zip( lengths,algs ) ]
    #stdev
    stdev   = std( ctgsCov )
    #max & min
    maxCov  = max( ctgsCov )
    maxCovC = contigs[ ctgsCov.index(maxCov) ]
    minCov  = min( ctgsCov )
    minCovC = contigs[ ctgsCov.index(minCov) ]

    #print info
    lines = "%s contigs/chromosomes\nReads per base:\t%.2f (+-%.2f)\nRange: %.2f [%s] - %.2f [%s]\n" % ( len(contigs),meanCov,stdev,minCov,minCovC,maxCov,maxCovC )
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
    coverage( o.fpath,o.outfn )
  
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
