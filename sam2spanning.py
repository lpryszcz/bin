#!/usr/bin/env python2
desc="""Parse BWA MEM sam output and report reads matching two contigs.
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 2/08/2013
"""

import argparse, os, sys
from datetime import datetime

def report_spanning(outfiles, sams, lenTh, verbose):
    """Take care only about the split_reads length
    and wether aligned twice.
M00724:10:000000000-A2EKA:1:1101:14333:1437     0       chrIV   233502  218     238M237S        *       0       0       TATCTAAAACATTTATAGGTTTGATTGTCATTCCTATTGTGGGTAATGCCGCAGAGCATGTCACTTCAGTCTTGGTGGCCATGAAGGATAAGATGGATCTGGCGCTAGGTGTTGCCATCGGTTCCTCTTTACAAGTTGCCTTATTTGTTACACCATTCATGGTTCTTGTGGGCTGGATGATCGCTGTTCCAATGACGCTAAATTTCTCACCTTTTGAAACCGCTACTCTTTTTATTGCGATGAGGTTGACGTAAAACAGCTTTCCACTGTGGCCACATTTCATCAACTTTCGAGGTAATTTCCAAATTCAGAAGTACAAGCCGCAAAGAATGCAGCACTCAGAGGAAAAGAAAAAGAAGAAGAGAAGGTCGGTTCTAGCCCACTAATACAGAAGTTGAAAAATGAAGATATCGAGTCTATCAAATGCAGGAACAATAATCTTTTAGACGGTAAAAAATTACTCGTGGAGGCAGAA     ???<?BBBDDDBDB<BFCFFCF?CCFFHFBFGHFF?CF/AE>E>E=0DDCCC5>CC=FHDCBGH=?DCGEAADBGFHEFDHH=AACDEFECFFHFFC.7EDG?)+)))?.@-=CDF-CFEED++,6=DDD=6=DDDEEC?=B,=BCC?C===BE4=?B,=BB??,?:A?*?A)08?A:A*?AA8*8.A?A**0:::A?8?8C**00*0*08:C:::8?**082;''0:*0*0:*0*0*0.8)0:*::?**).000:*00*0**0:*?:?A:*.8*08??##0.#*#:/#0?06*5.#0.06/(/6(6/(;E6.((-'.):8)@E@8@8*9*88;@DDD9EEEEDDD9EDDEED9ED@+EEED=CD<)74)+,EEDEA5DEEE@EEEEEE-CC5--55-C=C-FFC@A5+A5>@-C@DAA..FA8-A.FDC.9..A.DDA@9C.+>@+@A8,CBA8//CC@=<+-5-5-5/,,=,<     AS:i:226        XS:i:0  XF:i:3  XE:i:5  NM:i:3
M00724:10:000000000-A2EKA:1:1101:14333:1437     16      chrXIV  168815  212     237M238S        *       0       0       TTCTGCCTCCACGAGTAATTTTTTACCGTCTAAAAGATTATTGTTCCTGCATTTGATAGACTCGATATCTTCATTTTTCAACTTCTGTATTAGTGGGCTAGAACCGACCTTCTCTTCTTCTTTTTCTTTTCCTCTGAGTGCTGCATTCTTTGCGGCTTGTACTTCTGAATTTGGAAATTACCTCGAAAGTTGATGAAATGTGGCCACAGTGGAAAGCTGTTTTACGTCAACCTCATCGCAATAAAAAGAGTAGCGGTTTCAAAAGGTGAGAAATTTAGCGTCATTGGAACAGCGATCATCCAGCCCACAAGAACCATGAATGGTGTAACAAATAAGGCAACTTGTAAAGAGGAACCGATGGCAACACCTAGCGCCAGATCCATCTTATCCTTCATGGCCACCAAGACTGAAGTGACATGCTCTGCGGCATTACCCACAATAGGAATGACAATCAAACCTATAAATGTTTTAGATA     <,=,,/5-5-5-+<=@CC//8ABC,8A@+@>+.C9@ADD.A..9.CDF.A-8AF..AAD@C-@>5A+5A@CFF-C=C-55--5CC-EEEEEE@EEED5AEDEE,+)47)<DC=DEEE+@DE9DEEDDE9DDDEEEE9DDD@;88*9*8@8@E@)8:).'-((.6E;(/6(6/(/60.0#.5*60?0#/:#*#.0##??80*8.*:A?:?*:0**0*00*:000.)**?::*:0)8.0*0*0*:0*0*:0'';280**?8:::C:80*0*00**C8?8?A:::0**A?A.8*8AA?*A:A?80)A?*?A:?,??BB=,B?=4EB===C?CCB=,B=?CEEDDD=6=DDD=6,++DEEFC-FDC=-@.?)))+)?GDE7.CFFHFFCEFEDCAA=HHDFEHFGBDAAEGCD?=HGBCDHF=CC>5CCCDD0=E>E>EA/FC?FFHGFBFHFFCC?FCFFCFB<BDBDDDBBB?<???     AS:i:197        XS:i:0  XF:i:3  XE:i:4  NM:i:10    
    """
    if len(sams) < 2:
        return

    qname, tname1, seq1, qual1 = sams[0]
    qname, tname2, seq2, qual2 = sams[0]
    
    #check split length
    if len(seq1) < lenTh or len(seq2) < lenTh:
        sys.stderr.write( " %s < %s bases: CIGAR (%s,%s)!\n" % (qname,lenTh,"".join(cigar1),"".join(cigar2)) )
        return
        
    #report seqs and qualities into seperate files
    outfiles[0].write(">%s\n%s\n" % (qname, seq1))
    outfiles[1].write(">%s\n%s\n" % (qname, qual1))

def sam2spanning(handle, outfiles, mapqTh, lenTh, verbose):
    """Convert split algs into paired-end.
    i.e
SRR942191.59	16	xfSc0002356	1	60	152S115M1D18M1D10M	*	0	0	NNNCCTACTCCCCTGGTGTCGCCGTTGGCAGTCTCAGGCGGTCTTAGAAATCGCCAGTTCGTTGGTACTTATACTTATTTCGCCCAAATATTAATGGCTCGTTGTTATCTTCTTCTGCAAGAGAGGAACCAAGGGTTATTTCCTATTTTGTACAGAGCTTAGGCTTTCCAAGTCATTGTCGTGATCATACAAAGTGCTTGATGGTGATTGTTCGGTTCCGGTACTGTTCAATGTGCTGTGCAGTTTGTCACGACCATTTGAGTTCAATTTGCAAGAGGTATAGGCTTCTTGTAAA	!!!115:38,,,,----22151-,--55::755779664///5482?884:6985721,.,--//,-38=A==<..5444-0///44:AAABABEGHIFFFFFFIIIIIIIIIIIIIIIIIIF:97@@BBBB@55227;??;;3D==>>HHIIIIFFEHII====9999IIGGHE@BBFIIIIIIIFFDFHHHIIIIIIIIFDDEEEDFBBBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIHHHIIIIIIIIII666IIIIIIIIHHHIIIIIIIIIIIIIIIIIIIII@@@@	NM:i:11	AS:i:88	XS:i:0	SA:Z:xfSc0000404,9077,-,82S54M159S,60,0;
SRR942191.59	2064	xfSc0000404	9077	60	82H54M159H	*	0	0	CCCAAATATTAATGGCTCGTTGTTATCTTCTTCTGCAAGAGAGGAACCAAGGGT	///44:AAABABEGHIFFFFFFIIIIIIIIIIIIIIIIIIF:97@@BBBB@552	NM:i:0	AS:i:54	XS:i:0	SA:Z:xfSc0002356,1,-,152S115M1D18M1D10M,60,11;    
    """
    pname = ""
    sams  = []
    for sam in handle:
        #skip header and empty lines
        if sam.startswith("@") or not sam.strip():
            continue
            
        #for reference look at http://samtools.sourceforge.net/samtools.shtml
        qname,flag,tname,pos,mapq,cigar,mrnm,mpos,tlen,seq,qual = sam.split("\t")[:11]
        #skip algs with low quality
        if int(mapq) < mapqTh:
            continue
            
        #process stored algs if new read
        if qname != pname:
            report_spanning(outfiles, sams, lenTh, verbose)
            pname = qname
            sams  = []

        #store sam fields
        sams.append((qname, tname, seq, qual))

    #process very last entry
    report_spanning(outfiles, sams, lenTh, verbose)

def main():

    usage  = "bwa mem -t 3 contigs.fa ../../_archives/SRR942191.fastq.gz | %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", "--input",   default=sys.stdin, type=argparse.FileType('r'),
                        help="input stream        [stdin]")
    parser.add_argument("-o", "--outbase", required=True, #default=sys.stdout, type=argparse.FileType('w'),
                        help="output file name    [mandatory]")
    parser.add_argument("-q", "--mapq",    default=10,  type=int,
                        help="min mapping quality [%(default)s]")
    parser.add_argument("-l", "--minLen",  default=50, type=int,
                        help="min fragment length [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #define outnames
    fn1 = "%s.fasta" % o.outbase
    fn2 = "%s.fasta.qual" % o.outbase
    #get opened files
    outfiles = [open(fn1,"w"),open(fn2,"w")]
    #
    sam2spanning(o.input, outfiles, o.mapq, o.minLen, o.verbose)
        
if __name__=='__main__': 
    t0 = datetime.now()
    main()
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)    