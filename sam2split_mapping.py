#!/usr/bin/env python
desc="""Parse BWASW sam output and convert split matches into pair-end fastq.
NOTE: sam has to be sorted by read name (default BWASW output)
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 18/1/2013
"""

import argparse, os, re, sys
from Bio      import Seq
from datetime import datetime

cigarPat = re.compile('\d+\w')

def process_algs( outfiles,sams,lenTh,verbose ):
    """Take care only about the split_reads length
    and wether aligned twice.
M00724:10:000000000-A2EKA:1:1101:14333:1437     0       chrIV   233502  218     238M237S        *       0       0       TATCTAAAACATTTATAGGTTTGATTGTCATTCCTATTGTGGGTAATGCCGCAGAGCATGTCACTTCAGTCTTGGTGGCCATGAAGGATAAGATGGATCTGGCGCTAGGTGTTGCCATCGGTTCCTCTTTACAAGTTGCCTTATTTGTTACACCATTCATGGTTCTTGTGGGCTGGATGATCGCTGTTCCAATGACGCTAAATTTCTCACCTTTTGAAACCGCTACTCTTTTTATTGCGATGAGGTTGACGTAAAACAGCTTTCCACTGTGGCCACATTTCATCAACTTTCGAGGTAATTTCCAAATTCAGAAGTACAAGCCGCAAAGAATGCAGCACTCAGAGGAAAAGAAAAAGAAGAAGAGAAGGTCGGTTCTAGCCCACTAATACAGAAGTTGAAAAATGAAGATATCGAGTCTATCAAATGCAGGAACAATAATCTTTTAGACGGTAAAAAATTACTCGTGGAGGCAGAA     ???<?BBBDDDBDB<BFCFFCF?CCFFHFBFGHFF?CF/AE>E>E=0DDCCC5>CC=FHDCBGH=?DCGEAADBGFHEFDHH=AACDEFECFFHFFC.7EDG?)+)))?.@-=CDF-CFEED++,6=DDD=6=DDDEEC?=B,=BCC?C===BE4=?B,=BB??,?:A?*?A)08?A:A*?AA8*8.A?A**0:::A?8?8C**00*0*08:C:::8?**082;''0:*0*0:*0*0*0.8)0:*::?**).000:*00*0**0:*?:?A:*.8*08??##0.#*#:/#0?06*5.#0.06/(/6(6/(;E6.((-'.):8)@E@8@8*9*88;@DDD9EEEEDDD9EDDEED9ED@+EEED=CD<)74)+,EEDEA5DEEE@EEEEEE-CC5--55-C=C-FFC@A5+A5>@-C@DAA..FA8-A.FDC.9..A.DDA@9C.+>@+@A8,CBA8//CC@=<+-5-5-5/,,=,<     AS:i:226        XS:i:0  XF:i:3  XE:i:5  NM:i:3
M00724:10:000000000-A2EKA:1:1101:14333:1437     16      chrXIV  168815  212     237M238S        *       0       0       TTCTGCCTCCACGAGTAATTTTTTACCGTCTAAAAGATTATTGTTCCTGCATTTGATAGACTCGATATCTTCATTTTTCAACTTCTGTATTAGTGGGCTAGAACCGACCTTCTCTTCTTCTTTTTCTTTTCCTCTGAGTGCTGCATTCTTTGCGGCTTGTACTTCTGAATTTGGAAATTACCTCGAAAGTTGATGAAATGTGGCCACAGTGGAAAGCTGTTTTACGTCAACCTCATCGCAATAAAAAGAGTAGCGGTTTCAAAAGGTGAGAAATTTAGCGTCATTGGAACAGCGATCATCCAGCCCACAAGAACCATGAATGGTGTAACAAATAAGGCAACTTGTAAAGAGGAACCGATGGCAACACCTAGCGCCAGATCCATCTTATCCTTCATGGCCACCAAGACTGAAGTGACATGCTCTGCGGCATTACCCACAATAGGAATGACAATCAAACCTATAAATGTTTTAGATA     <,=,,/5-5-5-+<=@CC//8ABC,8A@+@>+.C9@ADD.A..9.CDF.A-8AF..AAD@C-@>5A+5A@CFF-C=C-55--5CC-EEEEEE@EEED5AEDEE,+)47)<DC=DEEE+@DE9DEEDDE9DDDEEEE9DDD@;88*9*8@8@E@)8:).'-((.6E;(/6(6/(/60.0#.5*60?0#/:#*#.0##??80*8.*:A?:?*:0**0*00*:000.)**?::*:0)8.0*0*0*:0*0*:0'';280**?8:::C:80*0*00**C8?8?A:::0**A?A.8*8AA?*A:A?80)A?*?A:?,??BB=,B?=4EB===C?CCB=,B=?CEEDDD=6=DDD=6,++DEEFC-FDC=-@.?)))+)?GDE7.CFFHFFCEFEDCAA=HHDFEHFGBDAAEGCD?=HGBCDHF=CC>5CCCDD0=E>E>EA/FC?FFHGFBFHFFCC?FCFFCFB<BDBDDDBBB?<???     AS:i:197        XS:i:0  XF:i:3  XE:i:4  NM:i:10    
    """
    if len(sams) != 2:
        return

    #compare cigars
    #for qname,cigarList,seq,qual in sams:
    ##237M,238S -> 238S,237M
    #it can be 267M,161S vs 266S,162M
    qname,cigar1,seq,qual = sams[0]
    qname,cigar2,seq,qual = sams[1]
    
    if   cigar1[0][-1]=="S":
        #if softclip from both ends #check which is larger
        if cigar1[-1][-1]=="S" and int(cigar1[-1][:-1]) > int(cigar1[0][:-1]):
            #if 3' then use 3'
            i = int(cigar1[-1][:-1])
        else:
            i = int(cigar1[0][:-1])
    elif cigar1[-1][-1]=="S":
        i = int(cigar1[-1][:-1])
    else:
        sys.stderr.write( " WARNING: query=%s CIGAR (%s,%s) handling problem!\n" % (qname,"".join(cigar1),"".join(cigar2)) )
        return
    #check split length
    if i<lenTh or len(seq)-i<lenTh:
        sys.stderr.write( " %s < %s bases: CIGAR (%s,%s)!\n" % (qname,lenTh,"".join(cigar1),"".join(cigar2)) )
        return
        
    #get FF fastq and report 
    outfiles[0].write( "@%s\n%s\n+\n%s\n" % ( qname,seq[:i],qual[:i] ) )
    outfiles[1].write( "@%s\n%s\n+\n%s\n" % ( qname,seq[i:],qual[i:] ) )

def split_mapping2pairs( handle,outfiles,mapqTh,lenTh,verbose ):
    """Convert split algs into paired-end.
    i.e
    """
    pname = ""
    sams  = []
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

        #237M238S -> 237M,238S
        cigarList = cigarPat.findall(cigar)
        #if mapped to reverse strand 
        if int(flag) & 16:
            #get reverse complement of seq and reverse quals
            seq   = str( Seq.Seq(seq).reverse_complement() )
            qual  = qual[::-1]
            #237M,238S -> 238S,237M
            cigarList.reverse() 
            
        #process stored algs if new read
        if qname!=pname:
            process_algs( outfiles,sams,lenTh,verbose )
            pname = qname
            sams  = []

        #store sam fields
        sams.append( (qname,cigarList,seq,qual) )

    #process very last entry
    process_algs( outfiles,sams,lenTh,verbose )

def main():

    usage  = "samtools view BAM [region(s)] | %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="input",   default=sys.stdin, type=argparse.FileType('r'),
                        help="input stream        [stdin]")
    parser.add_argument("-o", dest="outbase", required=True, #default=sys.stdout, type=argparse.FileType('w'),
                        help="output file name    [mandatory]")
    parser.add_argument("-q", dest="mapq",    default=10,  type=int,
                        help="min mapping quality [%(default)s]")
    parser.add_argument("-l", dest="minLen",  default=50, type=int,
                        help="min fragment length [%(default)s]")
  
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