#!/usr/bin/env python
"""
Buffered stdin sorting.
Don't see a way for this to work:/

Author:
l.p.pryszcz+git@gmail.com

Barcelona, 23/05/2012
"""

import sys
from datetime import datetime
from optparse import OptionParser
    
def main():
    usage  = "usage: %prog [options] fqA1 fqA2 [fqB1 fqB2]"
    desc   = """Buffered sorting stdin lines."""
    epilog = ""
    parser = OptionParser( usage=usage,version="%prog 1.0",description=desc,epilog=epilog ) 

    parser.add_option("-b", dest="buffsize", default="100", type=int,
                      help="buffer size (lines)    [%default]")
    parser.add_option("-v", dest="verbose", default=False, action="store_true" )

    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nArgs: %s\n" % ( o,args ) )

    buff = []
    l    = True
    while l:
        l = sys.stdin.readline()
        buff.append( l )
        if 
        
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )


