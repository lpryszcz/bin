#!/usr/bin/env python
"""Save all proteins from given phylome.

Author:
l.p.pryszcz+git@gmail.com

Barcelona, 29/05/2012
"""

import os, sys
from optparse import OptionParser
from datetime import datetime
from phylome  import _getConnection,_get_spcode

def main():
    usage = "usage: %prog [options]" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True

    parser.add_option("-p", dest="phyid", default=0, type=int,
                      help="define phylome id                 [mandatory]")
    parser.add_option("-s", dest="split",  default=False, action="store_true",
                      help="split fasta for ids from every line")
    parser.add_option("-v", dest="verbose",  default=True, action="store_false")
    
    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "%s\n" % ( str(o), ) )

    if not o.phyid:
        parser.error( "Specify mandatory parameters!" )

    #connect
    sys.stderr.write("Connecting...\n")
    p = _getConnection()
    
    #open common output
    sys.stderr.write("Fetching proteomes...\n")    
    if not o.split:
        outfn = "phylome_%s.fasta" % o.phyid
        out   = open( outfn,"w" )
    for proteome,pdata in p.get_proteomes_in_phylome(120)['proteomes'].iteritems():
        spcode,ver = proteome.split(".")
        taxid = pdata['taxid']
        sys.stderr.write(" %s (%s)       \r" % (spcode,proteome) )        
        #open output for each line if requested
        if o.split:
            outfn = "phylome_%s.%s.fasta" % ( o.phyid,proteome )
            out   = open( outfn,"w" )
        for id,data in p.get_seqs_in_genome(taxid,ver,filter_isoforms=True).iteritems():
            out.write( ">%s\n%s\n" % (id,data['seq']) )
        
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
