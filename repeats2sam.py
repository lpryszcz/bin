#!/usr/bin/env python2
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
l.p.pryszcz+git@gmail.com

Barcelona, 19/10/2012
"""

import argparse, os, sys, time
import subprocess
from datetime import datetime
from math     import exp
from Bio      import SeqIO

def get_sam(r, reads, s1, s2, rlen, gapfrac, inverted):
    """Return SAM enries for repeat
Col	Field	Description
1	QNAME	Query (pair) NAME
2	FLAG	bitwise FLAG
3	RNAME	Reference sequence NAME
4	POS	1-based leftmost POSition/coordinate of clipped sequence
5	MAPQ	MAPping Quality (Phred-scaled)
6	CIAGR	extended CIGAR string
7	MRNM	Mate Reference sequence NaMe ('=' if same as RNAME)
8	MPOS	1-based Mate POSistion
9	ISIZE	Inferred insert SIZE
10	SEQ	query SEQuence on the same strand as the reference
11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
12	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE    
    """
    flag1  = 67  #paired,mapped properly,1st in pair
    flag2  = 131 #paired,mapped properly,2nd in pair
    if inverted:
        flag1 += 32 #mate is reverse
        flag2 += 16 #read is reverse
        s2    -= rlen - 1 #change end to start
    chrid  = r.id
    rname  = "%s:%s-%s_%s" % ( chrid,s1,s2,rlen )
    if rname in reads:
        return "",rname        
    qual   = "I" * rlen # ord("I")=73 -> 73-33=40 PHRED+33
    mapq   = 255 #uniq match
    cigar  = "%sM" % rlen
    isize  = s2 - s1 + rlen
    seq1   = r.seq[ s1-1:s1+rlen-1 ]
    #skip if >
    if seq1.count('N')*1.0/rlen > gapfrac:
        return "",rname
    seq2   = r.seq[ s2-1:s2+rlen-1 ]
    sam    = ""
    sam   += "%s\t%s\t%s\t%s\t%s\t%s\t=\t%s\t%s\t%s\t%s\n" % ( rname,flag1,
                                    chrid,s1,mapq,cigar,s2,isize,seq1,qual )
    sam   += "%s\t%s\t%s\t%s\t%s\t%s\t=\t%s\t%s\t%s\t%s\n" % ( rname,flag2,
                                    chrid,s2,mapq,cigar,s1,isize,seq2,qual )
    return sam,rname
    
def chr2repeats(r,reads,out,lenth,len2distth,offset,maxdist,gapfrac,inverted,verbose,bufsize=-1):
    """
    bufsize: 0 no buffer; 1 buffer one line; -1 set system default.
    """
    fine = overlapping = farapart = 0
    #save chr to file
    fn = "/tmp/chr2repeat.%s.%s.tmp" % (r.id, time.time())
    tmpout = open( fn,"w" ); tmpout.write( ">%s\n%s\n" % (r.id,r.seq[offset:offset+maxdist]) ); tmpout.close()
    #get repeat-match command
    if not inverted:
        args = [ "repeat-match", "-f", "-n %s" % lenth, fn ]
    else:
        args = [ "repeat-match", "-n %s" % lenth, fn ]
        
    if verbose > 1:
        sys.stderr.write( "  %s\n" %  " ".join(args) )
    #run subprocess
    proc = subprocess.Popen( args,bufsize=bufsize,stdout=subprocess.PIPE,stderr=subprocess.PIPE )
    #parse stdout
    invrepeat = 0
    reads_new = set()
    for l in proc.stdout:
        l = l.strip()
        if not l or l.startswith(('#','Long','Start1','Exhaustive')):
            continue
        #get repeat start1 start2 and repeat length
        try:
            s1,s2,rlen = l.split()
            #control only for inverted, because for direct repeats -f is used
            if inverted:
                if s2.endswith('r'):
                    #strip r from start
                    s2 = s2[:-1]
                    invrepeat = 1
                else:
                    invrepeat = 0
            s1, s2, rlen = int(s1)+offset, int(s2)+offset, int(rlen)
        except:
            if verbose:
                sys.stderr.write( "Error: Cannot parse line  %s\n" % str(l.split()) )
            continue
            
        #check if passing distance criteria
        if s2 - s1 < rlen:
            overlapping += 1
            continue
        if rlen < 80 and s2 - s1 > len2distth[rlen]:
            farapart     += 1
            continue
        
        #get sam
        fine += 1
        sam,rname = get_sam(r, reads, s1, s2, rlen, gapfrac, invrepeat)
        reads_new.add(rname)
        if sam:
            out.write(sam)

    #remove tmp file
    os.unlink(fn)
    if verbose > 1:
        sys.stderr.write( "  %s repeats. fine : overlapping : farapart = %s : %s : %s \n" % ( fine+overlapping+farapart,fine,overlapping,farapart ) )
        
    return fine,overlapping,farapart,reads_new

def get_max_distances( lenth,uplimit,verbose ):
    """Return max distances for repeats having given length.
    """
    len2distth = {}
    for l in range( lenth,uplimit ):
        #lenlimit = 133.7480609953 * exp( 0.2011797391*l )
        lenlimit = 10*l**2.0 #2*l**2.0
        len2distth[l] = lenlimit
        sys.stderr.write( " %3i: %10i bp\n" % ( l,lenlimit ) )
    return len2distth        
    
def repeats2sam(handle, out, lenth, maxdist, gapfrac, inverted, verbose):
    """Report repeats in all chromosomes from multifasta.
    """
    #calculate distance boudaries for repeats in range lenth - 100
    if verbose:
        sys.stderr.write( "Max distances between repeats of length:\n" )        
    len2distth = get_max_distances( lenth,80,verbose )
    #process chr
    if verbose:
        sys.stderr.write( "Searching for repeats in chromosomes...\n" )

    statslines = "#chr\tlength\trepeats\tfine\toverlapping\tfarapart\n"
    for r in SeqIO.parse( handle,"fasta" ):
        if verbose:
            sys.stderr.write( " %s         \r" % r.id )
        #keep track of repeated entries
        reads = set()
        #print repeats
        for offset in range(0,len(r.seq),maxdist/2):
            fine,overlapping,farapart,reads = chr2repeats(r, reads, out, lenth, len2distth, offset, maxdist, gapfrac, inverted, verbose)
            #store stats
            statslines += "%s\t%s\t%s\t%s\t%s\t%s\n" % ( r.id,offset,fine+overlapping+farapart,fine,overlapping,farapart )

    #print stats if requested
    if verbose:
        sys.stderr.write( statslines )
        
def main():
    usage   = "%(prog)s [options] -i $ref | samtools view -SbuT $ref - | samtools sort - outfn"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )

    parser.add_argument("-v", dest="verbose",  default=False, action="store_true")
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="input",    default=sys.stdin, type=file, 
                        help="genome fasta      [stdin]" )
    parser.add_argument("-o", dest="output",   default=sys.stdout, type=argparse.FileType('w'), 
                        help="output            [stdout]" )    
    parser.add_argument("-l", dest="lenth",    default=8, type=int, 
                        help="min repeat length [%(default)s]" )
    parser.add_argument("-m", dest="maxdist",  default=10**5, type=int, 
                        help="max distance      [%(default)s]" )
    parser.add_argument("-g", dest="gapfrac",  default=0.1, type=float, 
                        help="max gap fraction [%(default)s]" )
    parser.add_argument("--inverted",          default=False, action="store_true",
                        help="report also inverted repeats" )
    
    o = parser.parse_args()
    
    if o.verbose: 
        sys.stderr.write( "Options: %s\n" % str(o) )

    #create outdir if not stdout output
    if type( o.output ) is file:
        outdir = os.path.dirname( o.output.name )
        if outdir and not os.path.isdir( outdir ):
            os.makedirs( outdir )

    #report repeats
    repeats2sam(o.input, o.output, o.lenth, o.maxdist, o.gapfrac, o.inverted, o.verbose)

if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
    