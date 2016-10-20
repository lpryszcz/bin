#!/usr/bin/env python
"""Report matches stats from blast.

Author:
l.p.pryszcz+git@gmail.com

Barcelona, 16/05/2012
"""

import os, sys
from optparse import OptionParser
from datetime import datetime
from genome_annotation import get_contig2size,parse_blast

def main():

    usage  = "usage: %prog [options] [ > out.bed ]"
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
    parser.add_option("-q", dest="qcov",   default=0, type=float,
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

    #get sizes of queries and targets
    q2len = get_contig2size( o.query  )
    t2len = get_contig2size( o.target )
    #get significant matches
    matches = parse_blast( o.infile,q2len,t2len,o.evalue,o.qcov,o.tcov,o.verbose )

    #parse matches
    print "#Query\tTarget\tIndentity\tAlignment length\tMismatches\tGaps\tQuery start\tQuery end\tTarget start\tTarget end\tE-value\tScore\tQuery aligned [%]\tTarget aligned [%]\t"
    for qlocus,tlocus,identity,algLen,mismatches,gaps,qstart,qend,tstart,tend,e,score,qcov,tcov in matches:
        '''#chr start end name score strand
        name = "%s:%s-%s" % ( tlocus,tstart,tend )
        #get strand
        strand = "+"
        if qstart>qend:
            strand = "-"
            qstart,qend = qend,qstart
        #define bed
        bed  = "%s\t%s\t%s\t%s\t%s\t%s\n" % ( qlocus,qstart-1,qend,name,score,strand ) '''
        out = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.1f\t%.1f\n" % (qlocus,tlocus,identity,algLen,mismatches,gaps,qstart,qend,tstart,tend,e,score,qcov*100,tcov*100, )
        sys.stdout.write( out )

    
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )


