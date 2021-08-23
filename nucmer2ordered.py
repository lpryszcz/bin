#!/usr/bin/env python2
"""Order contigs based on nucmer output.

Author:
l.p.pryszcz+git@gmail.com

Barcelona, 14/06/2012
"""

import os, sys
from optparse import OptionParser#,OptionGroup
from datetime import datetime
from genome_annotation import genome2dict,nucmer2list,_get_formatted_seq

def main():
    usage  = "usage: %prog [options]"
    desc   = """Order contigs based on nucmer output."""
    epilog = """Make sure, coords file is sorted by reference (show-coords -r).
Monoploid number (-x) has to be specified correctly. For details look at: http://en.wikipedia.org/wiki/Ploidy"""
    parser = OptionParser( usage=usage,version="%prog 1.0",description=desc,epilog=epilog ) 

    parser.add_option("-o", dest="outfn",  default="out",
                      help="output base name       [%default]")
    parser.add_option("-i", dest="coords", default="",
                      help="coords file name       [mandatory]")
    parser.add_option("-q", dest="query",  default="",
                      help="query file name        [mandatory]")
    parser.add_option("-r", dest="ref",    default="",
                      help="reference file name    [mandatory]")
    parser.add_option("-c", dest="qOverlap", default=0.05, type=float,
                      help="fract of query aligned [%default]")
    parser.add_option("-n", dest="haploid", default=2, type=int,
                      help="haploid number         [%default]")
    parser.add_option("-x", dest="monoploid", default=2, type=int,
                      help="monoploid number       [%default]")
    parser.add_option("-v", dest="verbose", default=False, action="store_true" )
  
    ( o, fnames ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nFastQ files: %s\n" % ( o,fnames ) )
        
    # check input files
    for fn in [ o.coords,o.query,o.ref ]:
        if not fn:
            parser.error( "Provide input file!" )
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )

    #load query genome
    query2fasta = genome2dict( o.query )
    ref2fasta   = genome2dict( o.ref )
        
    # load nucmer
    matches = nucmer2list( o.coords )

    # sort hits by ref position
    sort_hits( matches,query2fasta,ref2fasta,o.outfn,o.qOverlap,o.haploid,o.monoploid,o.verbose )

def sort_hits( matches,query2fasta,ref2fasta,outbase,qOverlapTh,haploid,monoploid,verbose ):
    """Return sorted multifasta of monoploid genomes.
    Contigs are sorted based on reference alignment.
    """
    #get best query to reference pairs
    ##prepare nested dictionary
    q2r = {}
    for q in query2fasta:
        q2r[q] = {}
    ##all query to ref
    for r,rStart,rStop,q,qStart,qStop,identity in matches:
        #for rStart,rStop,q,qStart,qStop,identity in matches[r]:
        qAligned = abs( qStop-qStart )
        #define if forward or reverse alg
        fwd = rev = 0        
        if qStop<qStart:
            rev = qAligned
        else:
            fwd = qAligned
        #store alg info
        if r not in q2r[q]:
            q2r[q][r] = [qAligned,rStart,rStop,fwd,rev]
            continue
        q2r[q][r][0] += qAligned
        if rStart<q2r[q][r][1]:
            q2r[q][r][1] = rStart
        if rStop>q2r[q][r][2]:
            q2r[q][r][2] = rStop
        #add fwd,rev
        q2r[q][r][3] += fwd
        q2r[q][r][4] += rev
        
    ##get best match for each query
    q2rBest = {}
    for q in query2fasta: #qSorted:
        refs = sorted( q2r[q].iteritems(), key=lambda x: q2r[q][x[0]][0], reverse=True )
        #skip contigs without a match or with too small fraction aligned
        if not refs or refs[0][1][0] < qOverlapTh * len(query2fasta[q]):
            continue
        q2rBest[q] = refs[0]
        print q,refs

    qOut = open( "%s.sorted.fa" % outbase,"w" )
    for r,rStart,rStop,q,qStart,qStop,identity in matches:
        if q not in q2rBest:
            continue
        #check if current r is the best for given q
        if r != q2rBest[q][0]:
            continue
        #pop given q from dictionary so it's saved only once
        qAligned,rStart,rStop,fwd,rev = q2rBest.pop(q)[1]
        #save sequence
        if fwd>rev:
            qOut.write( ">%s\n%s\n" % ( q, _get_formatted_seq( query2fasta[q] ) ) )
        #or it's reverse complement if more reverse aligned
        else:
            qOut.write( ">%s|rev\n%s\n" % ( q, _get_formatted_seq( query2fasta[q].reverse_complement() ) ) )
    return
    
    #separate monoploids
    ##how many monoploid sets
    noFiles = monoploid/haploid
    #define matches on reference global start and stops for each query
    matchesGlobal = {}
    for q,data in q2rBest.iteritems():
        r = data[0]
        qAligned,rGlobStart,rGlobStop = data[1]
        if r not in matchesGlobal:
            matchesGlobal[r]=[]
        matchesGlobal[r].append( (rGlobStart,rGlobStop,q) )
    #sort
    rOut = open( "%s.ref.sorted.fa" % outbase,"w" )
    qOut = open( "%s.query.sorted.fa" % outbase,"w" )
    for r in matchesGlobal:
        matchesGlobal[r].sort()

    for r in sorted( matchesGlobal.keys() ):
        rOut.write( ">%s\n%s\n" % ( r, _get_formatted_seq( ref2fasta[r] ) ) )
        for rGlobStart,rGlobStop,q in matchesGlobal[r]:
            print r,rGlobStart,rGlobStop,q
            qOut.write( ">%s\n%s\n" % ( q, _get_formatted_seq( query2fasta[q] ) ) )

    rOut.close()
    qOut.close()
    return
    for r in matchesGlobal:
        for rGlobStart,rGlobStop,q in matchesGlobal[r]:
            matchesInRegion = filter( lambda x: x[0]<rGlobStart<x[1] or x[0]<rGlobStop<x[1], matchesGlobal[r] )

            print r,rGlobStart,rGlobStop,q
            print matchesInRegion
            
    
    #save to file
    sys.stderr.write( "Going to save chromosomes in %s files:\n" % noFiles )
    outfiles=[]
    for i in range( noFiles ):
        outfn = "%s.%s.fa" % (outbase,i+1)
        sys.stderr.write( "  %s\n" % outfn )
        outfiles.append( open( outfn,"w" ) )

    #
    
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )

