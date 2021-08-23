#!/usr/bin/env python2
desc="""Generate sorted query ids list for mummerplot.
Example:
show-coords -r $ref.r.delta > $ref.r.coords
coords2ordered.py -i $ref.r.coords -f query.fasta  > query.ordered.txt
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

26/02/2013 Barcelona
"""

import argparse, commands, os, subprocess, sys
from datetime import datetime
#from genome_annotation import get_contig2size
from genome_annotation import genome2dict,nucmer2list

def coords2ordered(handle, out, fasta, minmatch, verbose):
    """ UNFINISHED"""
    #get contigs sizes
    contig2size = get_contig2size(fasta)
    #parse coords file
    queries = []
    contig2matches = {}
    for line in handle:
        #    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [FRM]  [TAGS]
        ts,te,qs,qe,tl,ql,identity,qframe,tframe,t,q = line.split('\t')
        ql = int(ql)
        #store match
        if q not in contig2matches:
            contig2matches[q] = {}
        if t not in contig2matches[q]:
            contig2matches[q][t] = [ ql ]
        else:
            contig2matches[q][t].append(ql)
            if sum(contig2matches[q][t]) > minmatch:
                queries.append([])
        #store ordered targets list
        if t not in refs:
            refs.append(t)
                    
def sort_hits( matches,query2fasta,out,verbose ):
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
    for q in query2fasta: 
        refs = sorted( q2r[q].iteritems(), key=lambda x: q2r[q][x[0]][0], reverse=True )
        #skip contigs without a match or with too small fraction aligned
        if not refs: # or refs[0][1][0] < qOverlapTh * len(query2fasta[q]):
            continue
        q2rBest[q] = refs[0]
        #print q,refs

    #qOut = open( "%s.sorted.fa" % outbase,"w" )
    for r,rStart,rStop,q,qStart,qStop,identity in matches:
        if q not in q2rBest:
            continue
        #check if current r is the best for given q
        if r != q2rBest[q][0]:
            continue
        #pop given q from dictionary so it's saved only once
        qAligned,rStart,rStop,fwd,rev = q2rBest.pop(q)[1]
        #get orientation
        orient = "+"
        if fwd<rev:
            orient = "-"
        #save sequence
        out.write( "%s\t%s\t%s\n" % ( q,len(query2fasta[q]),'+' ) )
 
def main():

    usage   = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )

    parser.add_argument("-v", dest="verbose", default=False, action="store_true")
    parser.add_argument("-i", dest="coords",  default=sys.stdin, type=file,
                        help="coords file name [stdin]")
    parser.add_argument("-f", dest="fasta",   required=True, type=file,
                        help="query fasta file" )
    parser.add_argument("-o", dest="out",  default=sys.stdout, type=argparse.FileType("w"),
                        help="output base name [stdout]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #load query genome
    query2fasta = genome2dict( o.fasta.name )

    # load nucmer
    matches = nucmer2list( o.coords.name )#; print matches

    # sort hits by ref position
    sort_hits( matches,query2fasta,o.out,o.verbose )                    
                    
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!   \n")
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
