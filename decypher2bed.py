#!/usr/bin/env python
"""Convert decypher output into BED.

Decypher has to be run with following output format:
[output format] tab #ncbi
[field] querylocus, targetlocus, rank, status, score, significance, querystart, queryend, querylength, targetstart, targetend, targetlength, gaps, matches, similarity

Author:
l.p.pryszcz@gmail.com

Barcelona, 25/04/2012
"""

import os, sys
from optparse import OptionParser
from datetime import datetime
from genome_annotation import parse_decypher

def main():

    usage  = "usage: %prog [options] [ > out.bed ]"
    desc   = """Decypher has to be run with following output format:
[output format] tab #ncbi
[field] querylocus, targetlocus, rank, status, score, significance, querystart, queryend, querylength, targetstart, targetend, targetlength, gaps, matches, similarity
    """
    epilog = ""
    parser = OptionParser( usage=usage,version="%prog 1.0",description=desc,epilog=epilog ) 

    parser.add_option("-i", dest="infile",  default="",
                      help="decypher output")
    parser.add_option("-e", dest="evalue", default=1e-05, type=float,
                      help="E-value cut-off [%default]" )
    parser.add_option("-q", dest="qcov",   default=0, type=float,
                      help="query coverage [%default]")
    parser.add_option("-t", dest="tcov",   default=0, type=float,
                      help="target coverage [%default]")
    parser.add_option("-s", dest="tsplit", default=3, type=int,
                      help="split target name by '|' and print s postition [%default]")        
    parser.add_option("-v", dest="verbose",  default=False, action="store_true" )    
  
    ( o, fnames ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nArgs: %s\n" % ( o,fnames ) )

    #check files
    for fn in ( o.infile, ):
        if not fn:
            parser.error( "Provide input file!" )
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )

    #get significant matches
    matches = parse_decypher( o.infile,o.evalue,o.qcov,o.tcov,o.verbose )

    #parse matches
    for qlocus,qstart,qend,tlocus,tstart,tend,score,e in matches:
        #chr start end name score strand
        name = "%s:%s-%s" % ( tlocus.split('|')[o.tsplit],tstart,tend )
        #get strand
        strand = "+"
        if qstart>qend:
            strand = "-"
            qstart,qend = qend,qstart
        #define bed
        bed  = "%s\t%s\t%s\t%s\t%s\t%s\n" % ( qlocus,qstart-1,qend,name,score,strand ) 
        sys.stdout.write( bed )
    
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
