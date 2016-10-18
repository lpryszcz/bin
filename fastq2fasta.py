#!/usr/bin/env python
desc="""Convert FastQ to FASTA.
Optionally store quals in separate file.
"""
epilog="""
l.p.pryszcz@gmail.com
Barcelona, 10/05/2013
"""
 
import argparse, gzip, os, sys
from datetime import datetime

def report(read, output, minLen, qualityTh, offset):
    """ """
    id, seq, spacer, quals = read
    #cut @ Ns
    if 'N' in seq:
        iN = seq.index('N')
        if iN<minLen:
            return 0
        seq = seq[:iN]
  	
    #cut @ low quality
    if qualityTh:
        qualsL = []
        for i in range(len(seq)):
            q = quals[i]
            if ord(q)-offset < qualityTh:
                break
            qualsL.append(q)
        #skip if too short	
        if i<minLen:
            return 0
        #trim seq
        seq = seq[:i+1]

    #store fasta				
    output.write('>%s\n%s\n' % (id[1:], "\n".join(seq[i:i+100] for i in range(0,len(seq),100))))
    return len(seq)

def fastq2fasta(input, output, minLen, qualityTh, offset, bases, verbose):
    """
    """
    #parse fastq
    i = 0
    read = []
    for line in input:
        #skip empty lines
        line = line.strip()
        if not line:
            continue
        #store read info
        read.append(line)
        #report reads
        if len(read)==4:
            i += report(read, output, minLen, qualityTh, offset)
            read = []
            # process up to bases
            if i>bases:
                break

def main():

    usage   = "%(prog)s [options] -v " 
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )

    parser.add_argument("-v", action="store_true", dest="verbose", default=False)
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-i", "--input",  default=sys.stdin, 
                        help="input file [stdin]")
    parser.add_argument("-o", "--output", default=sys.stdout, 
                        help="output file [stdout]")
    parser.add_argument("-l", "--minLen", default=0, type=int,
                        help="skip reads shorter than [%(default)s]" )
    parser.add_argument("-q", "--qualityTh", default=0, type=int,
                        help="read is clipped @ first base having PHRED quality lower than [%(default)s]" )
    #parser.add_argument("-q", dest="quals", default="", help="output quals to file [%(default)s]")
    parser.add_argument("--offset", default=33, type=int, 
                        help="quality encoding; PHRED+33 (Sanger) or PHRED+64 (Illumina/Solexa) [%(default)s]")
    parser.add_argument("-b", "--bases", default=float('inf'), type=float,
                        help="process up to b bases [%(default)s]" )

    o = parser.parse_args()

    #fastq2fasta( o.inFn,o.outFn,o.quals,o.phred,o.verbose )
    fastq2fasta(o.input, o.output, o.minLen, o.qualityTh, o.offset, o.bases, o.verbose)

if __name__=='__main__': 
    t0=datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt=datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
