#!/usr/bin/env python
"""
Parse multi-fasta file and report sequences without overlap
with already reported. Starts from the longest.

Author:
l.p.pryszcz+git@gmail.com

Dublin, 29/08/2012
"""

import os, sys
import commands
from optparse import OptionParser,OptionGroup
from datetime import datetime
from genome_annotation import genome2dict,_get_formatted_seq,parse_blat

def run_blat( fn1,fn2,minIdentity,verbose=1 ):
    """Execute blat and return output fname."""
    outfn = "%s.%s.psl" % ( os.path.basename(fn2),os.path.basename(fn1) )
    if not os.path.isfile( outfn ):
        cmd = "blat -q=dna -t=dna -extendThroughN -noHead -dots=100 -tileSize=11 -fastMap -minIdentity=%s %s %s %s" % (minIdentity,fn2,fn1,outfn)
        if verbose:
            sys.stderr.write( " %s\n" % cmd )
        os.system(cmd) 
    return outfn

def overlapping( c,added,matches,overlap,verbose ):
    """Execute blat and return output fname."""
    #Q,T,identity,Qalg,misM,Qgapbases,Qs,Qe,Ts,Te,strand,qcoverage,tcoverage
    #select only c aligned to bigger contigs than c and having overlap > threshold
    cmatches = filter( lambda x: x[0] == c and x[1] in added and x[11]>overlap, matches )# 
    if cmatches:
        #print cmatches
        return True
    return False

def main():
    usage  = "usage: %prog [options]"
    desc   = """Parse multi-fasta file and report sequences without overlap
with already reported sequences. Starts from the longest."""
    epilog = ""
    parser = OptionParser( usage=usage,version="%prog 1.0",description=desc,epilog=epilog ) 

    parser.add_option("-i", dest="infile",  
                      help="multi-fasta file       [mandatory]")
    parser.add_option("-m", dest="minIdentity",  default=90, type=int,
                      help="min identity           [%default]")
    parser.add_option("-o", dest="overlap",  default=0.3, type=float,
                      help="max overlap allowed    [%default]")
    parser.add_option("-v", dest="verbose", default=False, action="store_true" )
    
    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nArgs: %s\n" % ( o,args ) )

    for fn in [ o.infile, ]:
        if not fn:
            parser.error( "Provide input file!" )
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )

    #load fastas
    fastas = genome2dict( o.infile )

    #contigs by descending length
    contigs = sorted( fastas.keys(),key=lambda x: len(fastas[x]), reverse=True )

    #report non-overlapping
    i = 0
    added,skipped = set(), set()
    ##remove outfile if exists
    outfn = o.infile + ".collapsed_o%s_i%s.fa" % ( o.overlap,o.minIdentity )
    if os.path.isfile( outfn ):
        os.unlink( outfn )
    ##execute blat vs itself
    pslfn = run_blat( o.infile,o.infile,o.minIdentity,o.verbose )
    matches = parse_blat( pslfn,o.verbose,header=0,skipSelfMatches=1 )
    ##add contigs without overlap
    for c in contigs:
        i += 1
        if o.verbose:
            sys.stderr.write( " %3s %20s [ %7.2f kb]\n" % (i,c,len(fastas[c])/1000.0) )
        #get fasta entry
        fasta = ">%s\n%s\n" % (c,_get_formatted_seq(fastas[c]))
        #save contig if first or if no overlapping already processed
        if not added or not overlapping( c,added,matches,o.overlap,o.verbose ):
            added.add( c )
            out = open(outfn,"a"); out.write( fasta ); out.close()
        else:
            skipped.add( c )

    sys.stderr.write( "Selected %s [ %7.2f kb] out of %s [ %7.2f kb] contigs.\n" % ( len(added),sum([len(fastas[c]) for c in added])/10.0**3,len(fastas),sum([len(fastas[c]) for c in fastas])/10.0**3) )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )