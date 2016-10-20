#!/usr/bin/env python
"""Parses multiple blast files (for the same target) and print table
combining best matches to all target sequences from all query file.

Author:
l.p.pryszcz+git@gmail.com

Barcelona, 24/05/2012
"""

import os, sys
from optparse import OptionParser
from datetime import datetime
from genome_annotation import get_contig2size,parse_blast

def main():

    usage  = "usage: %prog [options] blastout1 [blastout2 ... blastoutN]  [ > out ]"
    desc   = """Blast has to be run with -m8."""
    epilog = ""
    parser = OptionParser( usage=usage,version="%prog 1.0",description=desc,epilog=epilog ) 

    parser.add_option("-k", dest="target",  default="",
                      help="target fasta")
    parser.add_option("-e", dest="evalue", default=1e-05, type=float,
                      help="E-value cut-off [%default]" )
    parser.add_option("-q", dest="qcov",   default=0, type=float,
                      help="query coverage  [%default]")
    parser.add_option("-t", dest="tcov",   default=0, type=float,
                      help="target coverage [%default]")
    parser.add_option("-s", dest="fnsplit", default=True, action="store_false",
                      help="split fnames    [%default]")        
    parser.add_option("-v", dest="verbose",  default=False, action="store_true" )    
  
    ( o, fnames ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nArgs: %s\n" % ( o,fnames ) )

    #check files
    for fn in fnames + [ o.target, ]:
        if not fn:
            parser.error( "Provide input file!" )
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )

    #get sizes of targets
    t2len = get_contig2size( o.target )

    #dict to store matches and list of targets
    s2matches = []
    targets   = sorted( t2len.keys() ) 
    
    #process all files
    samples = []
    for fn in fnames:
        #define sample name
        s = fn
        #split by dot if requested
        if o.fnsplit:
            s = fn.split(".")[0]
        samples.append( s )
        
        #define empty matches
        smatch = []
        for i in range( len(targets) ):
            smatch.append( [] )
            
        #get sizes of queries
        q2len = {}#get_contig2size( fn )        
        #get significant matches
        matches = parse_blast( fn,q2len,t2len,o.evalue,o.qcov,o.tcov,o.verbose )

        #parse matches
        for qlocus,tlocus,identity,algLen,mismatches,gaps,qstart,qend,tstart,tend,e,score,qcov,tcov in matches:
            i = targets.index( tlocus )
            #add match info if not match for given target
            if not smatch[i]:
                smatch[i] = ( qlocus,e,score,identity,tcov )
            #or better match found
            elif score > smatch[i][2]:
                smatch[i] = ( qlocus,e,score,identity,tcov )
        #store matches
        s2matches.append( smatch )

    #write header
    header = "Target"
    for s in samples:
        header += "\t%s\t" % s
    print header
    print "\t" + "identity [%]\tcoverage [%]\t" * len(samples)
    #write data
    for i in range( len(targets) ):
        line = targets[i]
        for smatch in s2matches:
            if smatch[i]:
                qlocus,e,score,identity,tcov = smatch[i]
            else:
                identity=tcov=0
            line += "\t%6.2f\t%6.2f" % ( identity,tcov*100 )
        print line

    
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )


