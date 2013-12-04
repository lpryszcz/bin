#!/usr/bin/env python
desc="""Report exact repeats in chromosomes.
Apply distance criteria as a function of repeat length:
 #f(x) = 133.7480609953 * EXP( 0.2011797391 * x )
 f(x) = 2x**2.0
 you may optimise this!
 there is a hard limit of 

Require MUMMER3 (repeat-match):
 -f only forward (ff) matches
 otherwise ff and fr matches:
   173938     182585r        9
   173938     200246         9  

"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 20/10/2012
"""

import argparse, os, sys
import subprocess
from datetime import datetime
from math     import exp, fabs 
from Bio      import SeqIO

def get_fasta( r,s1,s2,rlen,fraglen,verbose ):
    """Return fasta enries for repeat and flanking sequence 
    (flnk bp up- and downstream) and sequence of junction.
    """
    fasta  = ""
    #get repeat info
    chrid  = r.id
    rseq   = r.seq[ s1-1:s1+rlen-1 ]
    rname  = "%s:%s-%s_%s_%s" % ( chrid,s1,s2,rlen,rseq )
    #get flnk
    flnk  = ( fraglen - rlen ) / 2
    #start 1bp upstream if up- and downstream flanks are uneven:P
    oneoff  = 0 
    if ( fraglen - rlen ) % 2:
        oneoff = 1    
    for move in ( 0,-1*fraglen/4,1*fraglen/4 ):
        if flnk-3<fabs(move):
            continue
        #get seq and flanking bases around 1st repeat
        seq1   = r.seq[ s1-1-flnk-oneoff-move:s1+rlen-1+flnk-move ]
        #get seq and flanking bases around 2nd repeat    
        seq2   = r.seq[ s2-1-flnk-oneoff-move:s2+rlen-1+flnk-move ]
        #get junction sequence --> repeat, upstream 1st repeat and
        #downstream 2nd repeat
        seqj   = r.seq[ s1-1-flnk-oneoff-move:s1+rlen-1 ] + r.seq[ s2-1+rlen:s2+rlen-1+flnk-move ]
        #check if all seqs are of correct length
        if len(seq1)!=fraglen or len(seq2)!=fraglen or len(seqj)!=fraglen:
            if verbose:
                sys.stderr.write( "   Warning: Wrong seq length (%s:%s:%s) for %s [move: %s]!\n" % (len(seq1),len(seq2),len(seqj),rname,move) )
            continue
        #store as fasta
        fasta += ">%s.r1\n%s\n" % ( rname,seq1 )
        fasta += ">%s.r2\n%s\n" % ( rname,seq2 )
        fasta += ">%s.junction\n%s\n" % ( rname,seqj )    
    return fasta
    
def chr2repeats( r,out,lenth,maxlen,flnk,len2distth,verbose,bufsize=65536 ):
    """
    bufsize: 0 no buffer; 1 buffer one line; -1 set system default.
    """
    fine=overlapping=farapart=toolong=0
    #save chr to file
    fn = "/tmp/chr2repeat.tmp"
    tmpout = open( fn,"w" ); tmpout.write( ">%s\n%s\n" % (r.id,r.seq) ); tmpout.close()
    #get repeat-match command
    args = [ "repeat-match", "-f", "-n %s" % lenth, fn ] 
    if verbose:
        sys.stderr.write( "  %s\n" %  " ".join(args) )
    #run subprocess
    proc = subprocess.Popen( args,bufsize=bufsize,stdout=subprocess.PIPE )
    #parse stdout
    for l in proc.stdout:
        l = l.strip()
        if not l or l.startswith(('#','Long','Start1','Exhaustive')):
            continue
        #get repeat start1 start2 and repeat length
        s1,s2,rlen = l.split()
        try:
            s1,s2,rlen = int(s1),int(s2),int(rlen)
        except:
            if verbose:
                sys.stderr.write( "Error: Cannot parse line  %s\n" % str(l.split()) )
            continue
            
        #check if passing distance criteria
        if rlen > maxlen:
            toolong += 1
            continue
        if s2 - s1 < rlen:
            overlapping += 1
            continue
        if s2 - s1 > len2distth[rlen]:
            farapart     += 1
            continue
        
        #get sam
        fasta = get_fasta( r,s1,s2,rlen,flnk,verbose )
        if fasta:
            fine += 1
            out.write( fasta )

    #remove tmp file
    os.unlink( fn )
    if verbose:
        sys.stderr.write( "  %s repeats. fine : overlapping : farapart : toolong = %s : %s : %s : %s \n" % ( fine+overlapping+farapart+toolong,fine,overlapping,farapart,toolong ) )
        
    return fine,overlapping,farapart,toolong

def get_max_distances( lenth,uplimit,verbose ):
    """Return max distances for repeats having given length.
    """
    len2distth = {}
    for l in range( lenth,uplimit+1 ):
        #lenlimit = 133.7480609953 * exp( 0.2011797391*l )
        lenlimit = 2*l**2.0
        len2distth[l] = lenlimit
        sys.stderr.write( " %3i: %10i bp\n" % ( l,lenlimit ) )
    return len2distth            
    
def repeats2fasta( handle,out,lenth,maxlen,flnk,verbose ):
    """Report repeats in all chromosomes from multifasta.
    """
    #calculate distance boudaries for repeats in range lenth - 100
    if verbose:
        sys.stderr.write( "Max distances between repeats of length:\n" )        
    len2distth = get_max_distances( lenth,maxlen,verbose )
    #process chr
    if verbose:
        sys.stderr.write( "Scanning chromosomes for repeats...\n" )

    statslines = "#chr\tlength\trepeats\tfine\toverlapping\tfarapart\ttoolong\n"
    for r in SeqIO.parse( handle,"fasta" ):
        if verbose:
            sys.stderr.write( " %s         \n" % r.id )
        #print repeats
        fine,overlapping,farapart,toolong = chr2repeats( r,out,lenth,maxlen,flnk,len2distth,verbose )
        #store stats
        statslines += "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ( r.id,len(r.seq),
                fine+overlapping+farapart+toolong,fine,overlapping,farapart,toolong )

    #print stats if requested
    if verbose:
        sys.stderr.write( statslines )
        
def main():
    usage   = "usage: %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )

    parser.add_argument("-v", dest="verbose",  default=False, action="store_true")
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="input",    default=sys.stdin, type=file, 
                        help="genome fasta      [stdin]" )
    parser.add_argument("-o", dest="output",   default=sys.stdout, type=file, #argparse.FileType('w'), 
                        help="output            [stdout]" )    
    parser.add_argument("-l", dest="lenth",    default=8, type=int, 
                        help="min repeat length [%(default)s]" )
    parser.add_argument("-m", dest="maxlen",   default=30, type=int, 
                        help="max repeat length [%(default)s]" )
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

    #check maxlen 
    if o.lenth  >  o.maxlen:
        parser.error( "Min repeat length (-l) has to be smaller than max repeat length (-m)!" )
    if o.maxlen >= o.flnk:
        parser.error( "Reported region (-f) has to be bigger than repeat!" )
        
    #report repeats
    repeats2fasta( o.input,o.output,o.lenth,o.maxlen,o.flnk,o.verbose )

if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )