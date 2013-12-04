#!/usr/bin/env python
"""Parses multiple files and find common rows between these.

Author:
l.p.pryszcz@gmail.com

Barcelona, 17/04/2012
"""

import os, sys
from optparse import OptionParser
from datetime import datetime

def load_gene2pfam( pfam,eTh=1e-05,sTh=0 ):
    """Return dictionary of PFAM matches for
    every gene. { gene: [ ( tname,tacc,desc ) ] , }
    """
    gene2pfam = {}
    # target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description
    for l in open( pfam ):
        l = l.strip()
        if not l or l.startswith("#"):
            continue
        tname,tacc,qname,qacc,e,score,bias,bde,bdscore,bdbias,exp,reg,clu,ov,env,dom,rep,inc = l.split()[:18]
        desc    = " ".join( l.split()[18:] )
        qname   = qname.split('|')[-1]
        # filter out based on E-value and score
        e,score = float(e),float(score)
        if e>eTh or score<sTh:
            continue
        #store
        if not qname in gene2pfam:
            gene2pfam[qname] = []
        matchTuple = ( tname,tacc,desc )
        gene2pfam[qname].append( matchTuple )
    return gene2pfam

def add_comment_to_gtf( gtf,gene2pfam ):
    """Insert comments to gtf entries."""
    import urllib
    id = None
    for line in open( gtf ):
        line=line.strip() 
        if line.startswith('#') or not line: 
            continue
      
        line_data = line.split('\t')
        if len(line_data)<8: 
            continue # skip incorrect lines
    
        contig,source,feature,start,end,score,strand,frame,comments=line_data
    
        description={}
        for atr_value in comments.split(';'):
            atr_value = atr_value.strip()
            if not atr_value:
                continue
            atr,value=atr_value.split()
            if atr=='transcript_id':
                id = value.strip('"')

        # add description to comments
        pfam=pfams=''
        if id in gene2pfam:
            for name,acc,desc in gene2pfam[id]:
                pfams += "%s [%s]; " % ( desc,acc )
            pfam  = ' pfam "%s"' % urllib.quote( pfams[:-2] )
        
        sys.stdout.write( "%s%s\n" % ( line,pfam ) )
        
def main():
    
    usage = "usage: %prog [options]" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True

    parser.add_option("-g", dest="gtf",
                      help="genome annotation gtf/gff" )
    parser.add_option("-i", dest="pfam",
                      help="pfam tblout file")
    parser.add_option("-e", dest="evalue", default=1e-05, type=float,
                      help="E-value cut-off [%default]")
    parser.add_option("-s", dest="score", default=0    , type=float,
                      help="score cut-off   [%default]")
    parser.add_option("-v", dest="verbose",  default=True, action="store_false")
    
    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "%s\n" % ( str(o), ) )

    # check if files exists
    for fn in ( o.gtf,o.pfam ):
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" )

    #get gene2pfam
    gene2pfam = load_gene2pfam( o.pfam,o.evalue,o.score )
    
    # process gtf
    add_comment_to_gtf( o.gtf,gene2pfam )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
