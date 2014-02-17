#!/usr/bin/env python
"""

Author:
l.p.pryszcz@gmail.com

Dublin, 2/07/2012
"""

import commands, gzip, os, sys, subprocess
from optparse import OptionParser,OptionGroup
from datetime import datetime
from Bio      import Seq

def main():
    usage  = "usage: %prog [options] "
    desc   = ""
    epilog = ""
    parser = OptionParser( usage=usage,version="%prog 1.0",description=desc,epilog=epilog ) 

    parser.add_option("-o", dest="outdir",  default="simple_assembly",
                      help="output directory       [%default]")
    parser.add_option("-v", dest="verbose", default=False, action="store_true" )

    '''bwto = OptionGroup(parser, "Aligner options")
    bwto.add_option("-r", dest="ref",      default="",
                      help="reference fasta        [mandatory]")
    bwto.add_option("-b", dest="bwtopts",  default="--very-fast-local",
                      help="""bowtie2 options        [%default]\nNote, -k 1 --mm is added automatically; NEVER (!) use -p as it mess up read order in output!\nfor more: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml""")
    parser.add_option_group(bwto)'''
  
    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nArgs: %s\n" % ( o,args ) )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )


