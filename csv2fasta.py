#!/usr/bin/env python2
"""Retrieve fasta sequences from multifasta for given ids loaded from input file.

Author:
l.p.pryszcz+git@gmail.com

Barcelona, 29/05/2012
"""

import os, sys
from optparse import OptionParser
from datetime import datetime
from Bio      import SeqIO
from genome_annotation import genome2dict

def main():
    usage = "usage: %prog [options]" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True

    parser.add_option("-i", dest="input",  default="",
                      help="input file with tab-separated ids [mandatory]" )
    parser.add_option("-f", dest="fasta",  default="",
                      help="multifasta file                   [mandatory]" )
    parser.add_option("-o", dest="out",   default="out/out",
                      help="output fname                      [%default]" )
    parser.add_option("-s", dest="split",  default=False, action="store_true",
                      help="split fasta for ids from every line")
    parser.add_option("-v", dest="verbose",  default=True, action="store_false")
    
    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "%s\n" % ( str(o), ) )

    for f in ( o.input,o.fasta ):
        if not f:
            parser.error( "Specify mandatory parameters!" )
        if not f.isdigit and not os.path.isfile( f ):
            parser.error( "No such file: %s" % f )

    #check if outdir exists
    outdir = os.path.dirname(o.out)
    if not os.path.isdir( outdir ):
        os.makedirs( outdir )

    #get id2fasta
    sys.stderr.write( "Loading multifasta...\n" )
    id2fasta = genome2dict( o.fasta )

    sys.stderr.write( "Saving fastas...\n" )    
    #load ids
    i = 1
    #open common output
    if not o.split:
        outfn = "%s_%5i.fasta" % (o.out,i)
        outfn = outfn.replace(" ","0")
        out   = open( outfn,"w" )
    for l in open( o.input ):
        #open output for each line if requested
        if o.split:
            outfn = "%s_%5i.fasta" % (o.out,i)
            outfn = outfn.replace(" ","0")
            out   = open( outfn,"w" )
        for id in l.split():
            sys.stderr.write( " %s   \r" % id )
            if   id      in id2fasta:
                seq = id2fasta[id]
            elif id+"_1" in id2fasta:
                seq = id2fasta[id+"_1"]
            else:
                sys.stderr.write( "  No fasta for: %s\n" % id )
                continue
            out.write( ">%s\n%s\n" % (id,seq) )
        i += 1
        
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
