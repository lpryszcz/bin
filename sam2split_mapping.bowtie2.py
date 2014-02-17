#!/usr/bin/env python
desc="""Parse bowtie2 local bam and report split reads for partially aligned reads.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 18/1/2013
"""

import argparse, os, re, sys
from Bio      import Seq
from datetime import datetime

cigarPat = re.compile('\d+\w')
patM = re.compile(r'\d+M')

def process_alg( outfiles,qname,cigar,seq,qual,lenTh,allowedCigs,verbose ):
    """Take care only about the split_reads length
    and whether aligned twice.
M00724:10:000000000-A2EKA:1:1101:14333:1437     0       chrIV   233502  218     238M237S        *       0       0       TATCTAAAACATTTATAGGTTTGATTGTCATTCCTATTGTGGGTAATGCCGCAGAGCATGTCACTTCAGTCTTGGTGGCCATGAAGGATAAGATGGATCTGGCGCTAGGTGTTGCCATCGGTTCCTCTTTACAAGTTGCCTTATTTGTTACACCATTCATGGTTCTTGTGGGCTGGATGATCGCTGTTCCAATGACGCTAAATTTCTCACCTTTTGAAACCGCTACTCTTTTTATTGCGATGAGGTTGACGTAAAACAGCTTTCCACTGTGGCCACATTTCATCAACTTTCGAGGTAATTTCCAAATTCAGAAGTACAAGCCGCAAAGAATGCAGCACTCAGAGGAAAAGAAAAAGAAGAAGAGAAGGTCGGTTCTAGCCCACTAATACAGAAGTTGAAAAATGAAGATATCGAGTCTATCAAATGCAGGAACAATAATCTTTTAGACGGTAAAAAATTACTCGTGGAGGCAGAA     ???<?BBBDDDBDB<BFCFFCF?CCFFHFBFGHFF?CF/AE>E>E=0DDCCC5>CC=FHDCBGH=?DCGEAADBGFHEFDHH=AACDEFECFFHFFC.7EDG?)+)))?.@-=CDF-CFEED++,6=DDD=6=DDDEEC?=B,=BCC?C===BE4=?B,=BB??,?:A?*?A)08?A:A*?AA8*8.A?A**0:::A?8?8C**00*0*08:C:::8?**082;''0:*0*0:*0*0*0.8)0:*::?**).000:*00*0**0:*?:?A:*.8*08??##0.#*#:/#0?06*5.#0.06/(/6(6/(;E6.((-'.):8)@E@8@8*9*88;@DDD9EEEEDDD9EDDEED9ED@+EEED=CD<)74)+,EEDEA5DEEE@EEEEEE-CC5--55-C=C-FFC@A5+A5>@-C@DAA..FA8-A.FDC.9..A.DDA@9C.+>@+@A8,CBA8//CC@=<+-5-5-5/,,=,<     AS:i:226        XS:i:0  XF:i:3  XE:i:5  NM:i:3
    """    
    if   cigar[0 ][-1] in allowedCigs and int(cigar[ 0][:-1])>=lenTh:
        i = int(cigar[ 0][:-1])
    elif cigar[-1][-1] in allowedCigs and int(cigar[-1][:-1])>=lenTh:
        i = int(cigar[-1][:-1])
    else:
        sys.stderr.write( " WARNING: query=%s CIGAR (%s) handling problem!\n" % (qname,"".join(cigar)) )
        return
    #check split length
    if i<lenTh or len(seq)-i<lenTh:
        sys.stderr.write( " %s < %s bases: CIGAR (%s)!\n" % (qname,lenTh,"".join(cigar)) )
        return
       
    #get FF fastq and report 
    outfiles[0].write( "@%s/1\n%s\n+\n%s\n" % ( qname,seq[:i],qual[:i] ) )
    outfiles[1].write( "@%s/2\n%s\n+\n%s\n" % ( qname,seq[i:],qual[i:] ) )

def split_mapping2pairs( handle,outfiles,mapqTh,lenTh,verbose ):
    """Convert split algs into paired-end.
    i.e
    """
    allowedCigs=set(["M","S"])
    for sam in handle:
        #directly print header lines(^@) 
        if sam.startswith("@"):
            continue
        #skip empty lines
        elif not sam.strip():
            continue
            
        #for reference look at http://samtools.sourceforge.net/samtools.shtml
        qname,flag,tname,pos,mapq,cigar,mrnm,mpos,tlen,seq,qual = sam.split("\t")[:11]
        #skip algs with low quality
        if int(mapq)<mapqTh:
            continue

        #get number of substitutions
        mtchs  = [ int(m[:-1]) for m in patM.findall( cigar ) ]
        #check at least 50bp is matched and 50bp is unmatched
        if not lenTh < sum(mtchs) < len(seq)-lenTh:
            continue
        #print sam
        #continue
        #237M238S -> 237M,238S
        cigarList = cigarPat.findall(cigar)#; print cigarList, seq
        #if mapped to reverse strand 
        if int(flag) & 16:
            #get reverse complement of seq and reverse quals
            seq   = str(Seq.Seq(seq).reverse_complement())
            qual  = qual[::-1]
            #237M,238S -> 238S,237M
            #cigarList.reverse()#; print 'reverse', cigarList, seq
        #mark as first or second from pair
        if   int(flag) & 64:
            qname += ".1"
        elif int(flag) & 128:
            qname += ".2"
        #process
        process_alg( outfiles,qname,cigarList,seq,qual,lenTh,allowedCigs,verbose )

def main():

    usage  = "samtools view BAM [region(s)] | %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="input",   default=sys.stdin, type=argparse.FileType('r'),
                        help="input stream        [stdin]")
    parser.add_argument("-o", dest="outbase", required=True, #type=argparse.FileType('w'),
                        help="output file name    [mandatory]")
    parser.add_argument("-q", dest="mapq",    default=10,  type=int,
                        help="min mapping quality [%(default)s]")
    parser.add_argument("-l", dest="minLen",  default=50, type=int,
                        help="min read length     [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #define outnames
    fn1 = "%s_1.fastq" % o.outbase
    fn2 = "%s_2.fastq" % o.outbase
    #get opened files
    outfiles = [ open(fn1,"w"),open(fn2,"w") ]
    #
    split_mapping2pairs( o.input,outfiles,o.mapq,o.minLen,o.verbose )
        
if __name__=='__main__': 
    t0 = datetime.now()
    main()
    dt = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )    