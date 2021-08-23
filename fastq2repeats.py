#!/usr/bin/env python2
desc="""Report exact matches in fastq reads for repeats.

Repeats are given as fasta.
Use repeats2fasta.py to generate it. 

"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 20/10/2012
"""

import argparse, gzip, os, sys
import subprocess
from datetime import datetime
from math     import exp
from Bio      import SeqIO

def repeats2dict( handle ):
    """Return dictionary of repeats.
    """
    seq2repeat={}
    repeat2cnt={}
    #open as gzip if endswith .gz
    if type(handle) is file and handle.name.endswith(".gz"):
        handle = gzip.open( handle.name )
    #parse lines
    for r in SeqIO.parse(handle,'fasta'):
        coord,rlen,rseq = r.id.split(".")[0].split("_")
        #you can generate some stats about repeats
        seq   = str(r.seq)
        seqrc = str(r.seq.reverse_complement())
        #store uniq seq or seqrc
        if seq not in seq2repeat:
            seq2repeat[seq  ] = []
            seq2repeat[seqrc] = []

        #add entry
        seq2repeat[seq  ].append( r.id )
        seq2repeat[seqrc].append( r.id )

        #add repeat2count entry
        repeat2cnt[r.id]=[0,0,0] # uniq,norm,all
    
    return seq2repeat,repeat2cnt

def seq2repeat( rseq2rname,repeat2cnt,seq,flnk,matches,verbose ):
    """
    """
    # repeat2cnt[r.id]=[0,0,0] # uniq,norm,all
    matchList = []
    for j in range( 0,len(seq)-flnk ):
        seqj = seq[j:j+flnk]
        if seqj in rseq2rname:
            matches += 1
            for rname in rseq2rname[seqj]:
                matchList.append( rname )
                repeat2cnt[rname][2] += 1
                
    #update uniq
    if len(matchList) == 1:
        repeat2cnt[matchList[0]][0] += 1

    #update norm
    for rname in matchList:
        repeat2cnt[rname][1] += 1.0 / len(matchList)
        
    return matchList,matches,repeat2cnt
    
def fastq2repeats( out,handle,repeats,flnk,readsfmt,verbose ):
    """
    """
    #load repeats
    if verbose:
        sys.stderr.write( "Loading repeats...\n" )
    rseq2rname,repeat2cnt = repeats2dict( repeats )
    if verbose:
        sys.stderr.write( " %s repeats loaded.\n" % len(repeat2cnt) )    
    #open as gzip if endswith .gz
    if type(handle) is file:
        if handle.name.endswith(".gz"):
            handle = gzip.open( handle.name )
            
    #process reads
    if verbose:
        sys.stderr.write( "Scanning reads...\n" )
        
    pI=umatches=matches=0
    step=10**5
    for i,r in enumerate( SeqIO.parse( handle,readsfmt ) ):
        #verbose
        if verbose and i>pI:
            pI += step
            sys.stderr.write( " %9i [ %9i : %9i ] Last read seen: %s    \r" % ( i,matches,umatches,r.id ) )
        #get seq and seqrc
        matchList,matches,repeat2cnt = seq2repeat( rseq2rname,repeat2cnt,str(r.seq),flnk,matches,verbose )
        #update stats        
        if matchList:
            umatches += 1
    #just to have an idea about number of reads
    if verbose:
        sys.stderr.write( " %s reads processed [ %s : %s ].\n" % (i,matches,umatches) )

    #save summary - readname - data for junction - data for repeat1 - data for repeat2
    lines  = "#\tjunction\t\trepeat1\t\trepeat2\t\t\n"
    lines += "#repeat name\tuniq\tnormalised\tall\n"
    for rname in sorted(repeat2cnt):
        uniq,norm,allm = repeat2cnt[rname]
        if rname.endswith( ".junction" ):
            lines += "%s\t" % rname[:-len(".junction")]
        lines += "%s\t%s\t%s\t" % ( uniq,norm,allm )
        if rname.endswith( ".r2" ):
            lines += "\n"
    out.write( lines )
    
def main():
    usage   = "usage: %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )

    parser.add_argument("-v", dest="verbose",  default=False, action="store_true")
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-r", dest="repeats",  default=None, type=file, 
                        help="repeat fasta      [mandatory]" )
    parser.add_argument("-i", dest="reads",    default=sys.stdin, type=file,
                        help="fasta/q file      [stdin]" )
    parser.add_argument("-t", dest="readsfmt", default="fastq", 
                        help="reads format      [%(default)s]" )    
    parser.add_argument("-o", dest="output",   default=sys.stdout, type=file, #argparse.FileType('w'), 
                        help="output            [stdout]" )    
    parser.add_argument("-f", dest="flnk",     default=40, type=int, 
                        help="report flanking   [%(default)s bp]" )
    
    o = parser.parse_args()
    
    if o.verbose: 
        sys.stderr.write( "Options: %s\n" % str(o) )

    #create outdir if not stdout output
    if type( o.output ) is file:
        outdir = os.path.dirname( o.output.name )
        if outdir and not os.path.isdir( outdir ):
            os.makedirs( outdir )

    #report repeats
    fastq2repeats( o.output,o.reads,o.repeats,o.flnk,o.readsfmt,o.verbose )

if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )