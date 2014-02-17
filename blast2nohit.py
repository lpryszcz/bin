#!/usr/bin/env python
"""Report queries with no match from blast.

Author:
l.p.pryszcz@gmail.com

Dublin, 17/07/2012
"""

import os, sys
from optparse import OptionParser
from datetime import datetime
from genome_annotation import get_contig2size, genome2dict, parse_blast, _get_formatted_seq
from Bio      import SeqIO

def main():

    usage  = "usage: %prog [options] [ 1> matches.table.txt ]"
    desc   = """Blast has to be run with -m8."""
    epilog = ""
    parser = OptionParser( usage=usage,version="%prog 1.0",description=desc,epilog=epilog ) 

    parser.add_option("-i", dest="infile",  default="",
                      help="blast output")
    parser.add_option("-j", dest="query",  default="",
                      help="query fasta")
    parser.add_option("-k", dest="target",  default="",
                      help="target fasta")
    parser.add_option("-e", dest="evalue", default=1e-05, type=float,
                      help="E-value cut-off [%default]" )
    parser.add_option("-q", dest="qcov",   default=0.3, type=float,
                      help="query coverage  [%default]")
    parser.add_option("-t", dest="tcov",   default=0, type=float,
                      help="target coverage [%default]")
    #parser.add_option("-s", dest="tsplit", default=3, type=int,
    #                  help="split target name by '|' and print s postition [%default]")        
    parser.add_option("-v", dest="verbose",  default=False, action="store_true" )    
  
    ( o, fnames ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nArgs: %s\n" % ( o,fnames ) )

    #check files
    for fn in ( o.infile,o.query,o.target ):
        if not fn:
            parser.error( "Provide input file!" )
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )

    #queries = get_
    #get sizes of queries and targets
    q2len = get_contig2size( o.query  )
    t2len = get_contig2size( o.target )
    #get significant matches
    matches = parse_blast( o.infile,q2len,t2len,o.evalue,0,0,o.verbose )

    #parse matches
    matches_collapsed = {}
    print "#Query\tTarget\tIndentity\tAlignment length\tMismatches\tGaps\tQuery start\tQuery end\tTarget start\tTarget end\tE-value\tScore\tQuery aligned [%]\tTarget aligned [%]\t"
    for qlocus,tlocus,identity,algLen,mismatches,gaps,qstart,qend,tstart,tend,e,score,qcov,tcov in matches:
        #add qlocus to matches
        if qlocus not in matches_collapsed:
            matches_collapsed[qlocus]={}
        if tlocus not in matches_collapsed[qlocus]:
            matches_collapsed[qlocus][tlocus]=[]
        #store data
        matches_collapsed[qlocus][tlocus].append( (algLen,mismatches,gaps,qstart,qend,tstart,tend,e,score,qcov,tcov) )

    #
    matched_queries = set()
    for qlocus in sorted( matches_collapsed.keys() ):
        for tlocus in sorted( matches_collapsed[qlocus].keys() ):
            qCov=tCov=0
            for algLen,mismatches,gaps,qstart,qend,tstart,tend,e,score,qcov,tcov in matches_collapsed[qlocus][tlocus]:
                qCov += qcov
                tCov += tcov

            if qCov<o.qcov or tCov<o.tcov:
                continue
            out = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.1f\t%.1f\n" % (qlocus,tlocus,identity,algLen,mismatches,gaps,qstart,qend,tstart,tend,e,score,qCov*100,tCov*100, )
            sys.stdout.write( out )
            matched_queries.add( qlocus )

    #get with no valid match
    sys.stderr.write( "Queries without valid matches:\n" )
    i = 0
    out = open( o.query + ".nomatch.fa","w" )
    for r in SeqIO.parse( open(o.query),"fasta" ):
        if r.id in matched_queries:
            continue
        i+=1
        line = "%s\t%s\t%s" % (i,r.id,len(r.seq))
        if r.id in matches_collapsed:
            line += "\t%s" % str( matches_collapsed[r.id] )
        sys.stderr.write( line+"\n" )
        #save fasta
        out.write( ">%s\n%s\n" % ( r.id,_get_formatted_seq( r.seq ) ) ) 
        
            
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )


