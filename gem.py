#!/usr/bin/env python2
"""Run gem-mapper

"""

import os, sys
from optparse import OptionParser
from datetime import datetime
import subprocess


def main():

    usage  = "usage: %prog [options] _archives/7A1_read?.fastq" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True

    parser.add_option("-d", dest="workingDir", default="",
                      help="change dir at first [%default]")
    parser.add_option("-I", dest="index",      default="../gem_index/F.oxysporum.fa",
                      help="gem index           [%default]")
    parser.add_option("-m", dest="mismatches", default=2, 
                      help="mismatches [%default]")
    parser.add_option("-v", dest="verbose",  default=False, action="store_true" )
  
    ( o, fnames ) = parser.parse_args()
    if o.verbose:
        print "Options: %s\nFiles: %s" % ( o,fnames )

    if o.workingDir:
        if o.verbose:
            print "Chanding directory: %s" % o.workingDir
        os.chdir( o.workingDir )
        
    if not fnames:
        parser.error( "Provide at least one fastq file!" )
    for fn in fnames:
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )

    cmd = "gem-mapper -q phred -I ../gem_index/F.oxysporum.fa -i _archives/7A1_read1_fastq -t1 -m5 -o gem/7A1.1 > gem/7A1.1.log 2>@1"
    
    
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
