#!/usr/bin/env python
desc="""Return sequences diverged by given percentage compares to input.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow, 4/03/2014
"""

import os, random, sys
from datetime import datetime
from Bio import SeqIO, Seq

aminos = 'ACDEFGHIKLMNPQRSTVWY'
nucleotides = 'ACGT'

def seq2diverged(seq, divergence, loh, aminos):
    """Return diverged sequence"""
    seqlist = list(seq)
    #number of position to change
    k = int(round(divergence*len(seqlist)))
    positions = random.sample(xrange(len(seqlist)), k)
    if loh:
        #what is the distribution of LOH sizes?? poisson?
        pass
    #change positions
    for i in positions:
        nb = seqlist[i]
        while nb == seqlist[i]:
            nb = random.choice(aminos)
        seqlist[i] = nb
    return "".join(seqlist)

def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-i", "--input",     default=sys.stdin, type=file, 
                        help="fasta file(s)   [stdin]")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-d", "--divergence",    default=0.01, type=float, 
                        help="divergence      [%(default)s]")
    parser.add_argument("--dna",             default=False, action='store_true',
                        help="DNA alphabet    [amino acids]")
    parser.add_argument("--loh",             default=0.00, type=float, 
                        help="level of LOH    [%(default)s]")
                        
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    for i, r in enumerate(SeqIO.parse(o.input, 'fasta'), 1):
        if o.verbose:
            sys.stderr.write(' %s\r'%i)
        if o.dna:
            alphabet = nucleotides
        else:
            alphabet = aminos
        seq = seq2diverged(r.seq, o.divergence, o.loh, alphabet)
        r.seq = Seq.Seq(seq)
        o.output.write(r.format('fasta'))

if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
