#!/usr/bin/env python
"""Reports GC of every chromosome/contig.

Author:
l.p.pryszcz@gmail.com

Barcelona, 15/06/2012
"""

import os, sys
from optparse import OptionParser#,OptionGroup
from datetime import datetime
from Bio      import SeqIO
from Bio.SeqUtils import GC
from genome_annotation import genome2dict

def main():
    usage  = "usage: %prog [options]"
    desc   = """Reports GC of every chromosome/contig."""
    epilog = ""
    parser = OptionParser( usage=usage,version="%prog 1.0",description=desc,epilog=epilog ) 

    #parser.add_option("-o", dest="outfn",  default="out",
    #                  help="output base name       [%default]")
    parser.add_option("-i", dest="infn", default="",
                      help="genome file name       [mandatory]")
    parser.add_option("-v", dest="verbose", default=False, action="store_true" )
  
    ( o, fnames ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nArgs: %s\n" % ( o,fnames ) )
        
    # check input files
    for fn in [ o.infn, ]:
        if not fn:
            parser.error( "Provide input file!" )
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )

    for r in SeqIO.parse( open(o.infn),"fasta" ):
        sys.stdout.write( "%s\t%s\t%s\n" % ( r.id,len(r.seq),GC(r.seq) ) )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
  
            