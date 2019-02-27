#!/usr/bin/env python
desc="""Return sequences diverged by given percentage compares to input.
LOH sizes are drawn from log-normal distribution. 
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Mizerow, 4/03/2014
"""

import os, random, sys
from datetime import datetime
from Bio import SeqIO, Seq
import numpy as np

aminos = 'ACDEFGHIKLMNPQRSTVWY'
nucleotides = 'ACGT'

JC69 = """ A C G T
A 1 1 1 1
C 1 1 1 1
G 1 1 1 1
T 1 1 1 1"""

# http://www.uky.edu/Classes/BIO/520/BIO520WWW/blosum62.htm
blossum62 = """ C	S	T	P	A	G	N	D	E	Q	H	R	K	M	I	L	V	F	Y	W
C	9	-1	-1	-3	0	-3	-3	-3	-4	-3	-3	-3	-3	-1	-1	-1	-1	-2	-2	-2
S	-1	4	1	-1	1	0	1	0	0	0	-1	-1	0	-1	-2	-2	-2	-2	-2	-3
T	-1	1	4	1	-1	1	0	1	0	0	0	-1	0	-1	-2	-2	-2	-2	-2	-3
P	-3	-1	1	7	-1	-2	-1	-1	-1	-1	-2	-2	-1	-2	-3	-3	-2	-4	-3	-4
A	0	1	-1	-1	4	0	-1	-2	-1	-1	-2	-1	-1	-1	-1	-1	-2	-2	-2	-3
G	-3	0	1	-2	0	6	-2	-1	-2	-2	-2	-2	-2	-3	-4	-4	0	-3	-3	-2
N	-3	1	0	-2	-2	0	6	1	0	0	-1	0	0	-2	-3	-3	-3	-3	-2	-4
D	-3	0	1	-1	-2	-1	1	6	2	0	-1	-2	-1	-3	-3	-4	-3	-3	-3	-4
E	-4	0	0	-1	-1	-2	0	2	5	2	0	0	1	-2	-3	-3	-3	-3	-2	-3
Q	-3	0	0	-1	-1	-2	0	0	2	5	0	1	1	0	-3	-2	-2	-3	-1	-2
H	-3	-1	0	-2	-2	-2	1	1	0	0	8	0	-1	-2	-3	-3	-2	-1	2	-2
R	-3	-1	-1	-2	-1	-2	0	-2	0	1	0	5	2	-1	-3	-2	-3	-3	-2	-3
K	-3	0	0	-1	-1	-2	0	-1	1	1	-1	2	5	-1	-3	-2	-3	-3	-2	-3
M	-1	-1	-1	-2	-1	-3	-2	-3	-2	0	-2	-1	-1	5	1	2	-2	0	-1	-1
I	-1	-2	-2	-3	-1	-4	-3	-3	-3	-3	-3	-3	-3	1	4	2	1	0	-1	-3
L	-1	-2	-2	-3	-1	-4	-3	-4	-3	-2	-3	-2	-2	2	2	4	3	0	-1	-2
V	-1	-2	-2	-2	0	-3	-3	-3	-2	-2	-3	-3	-2	1	3	1	4	-1	-1	-3
F	-2	-2	-2	-4	-2	-3	-3	-3	-3	-3	-1	-3	-3	0	0	0	-1	6	3	1
Y	-2	-2	-2	-3	-2	-3	-2	-3	-2	-1	2	-2	-2	-1	-1	-1	-1	3	7	2
W	-2	-3	-3	-4	-3	-2	-4	-4	-3	-2	-2	-3	-3	-1	-3	-2	-3	1	2	11"""

def load_sim_matrix(matrix=blossum62):
    """Return dict of weighted seq for substitution for each base"""
    lines = [l.strip().split() for l in matrix.split('\n') if l.strip()]
    base2subs = {}
    bases = lines[0]
    for i, l in enumerate(lines[1:]):
        b, scores = l[0], np.array(map(int, l[1:]))
        scores += -min(scores)+1
        scores[i] = 0
        base2subs[b] = "".join(_b*s for _b, s in zip(bases, scores))
    return base2subs

def get_heterozygous(positions, divergence):
    """Return fraction of chromosome being heterozygous"""
    #count hetero SNPs within 100bp
    hetero = []
    pp = 0.0
    for p in sorted(positions):
        if   p - pp <=100:
            if not hetero or len(hetero[-1])==2:
                hetero.append([pp,])
            pp = p
        elif hetero and len(hetero[-1])==1:
            hetero[-1].append(pp)
        pp = p
    if hetero and len(hetero[-1])==1:
        hetero[-1].append(pp)
    return sum(e-s for s, e in hetero)

def seq2diverged(seq, divergence, loh, lohSizes, base2subs, verbose):
    """Return diverged sequence"""
    seqlist = list(seq)
    #number of position to change
    k = int(round(divergence*len(seqlist)))
    positions = random.sample(xrange(len(seqlist)), k)
    #inclue LOHs
    if loh:
        lohs = []
        while get_heterozygous(positions, divergence)/len(seq) > 1-loh:
            # get LOH random start
            s = random.randint(0, len(seq)-1)
            # and random LOH length from negative binomial distribution
            lsize = 0
            # LOH has to be at least 2x larger than average distance between SNPs
            while lsize <= 2.0 / divergence:
                lsize = int(round(random.sample(lohSizes, 1)[0]))
            e = s + lsize
            lohs.append(lsize)
            # filter snp posiitons
            positions = filter(lambda p: p<s or p>e, positions)
        if verbose:
            hetero = 100.0 * get_heterozygous(positions, divergence)/len(seq)
            sys.stderr.write(" %.1f%s heterozygous with %s LOHs (%s - %sbp)\n" % \
                             (hetero, '%', len(lohs), min(lohs), max(lohs)))
    #change positions
    for i in positions:
        if seqlist[i].upper() in base2subs:
            seqlist[i] = random.choice(base2subs[seqlist[i].upper()])
    return "".join(seqlist)

def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
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
    parser.add_argument("--learn",    default=False, type=file, 
                        help="BED file to learn LOH sizes [%(default)s]")
    parser.add_argument("--power_a",         default=0.3, type=float, 
                        help="LOH power distribution a [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
    # load sim matrix
    if o.dna:
        base2subs = load_sim_matrix(matrix=JC69)
    else:
        base2subs = load_sim_matrix(matrix=blossum62)
    # get LOH sizes distibution
    lohSizes = []    
    if o.loh:
        lohSizes = np.random.lognormal(0, 1.2, 1e5) * 1e3
    # parse sequences
    if o.verbose:
        sys.stderr.write("Processing FastA...\n")
    for i, r in enumerate(SeqIO.parse(o.input, 'fasta'), 1):
        if o.verbose:
            sys.stderr.write('%s %s %sbp %s\n'%(i, r.id, len(r), " "*20))
        seq = seq2diverged(r.seq, o.divergence, o.loh, lohSizes, base2subs, o.verbose)
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
