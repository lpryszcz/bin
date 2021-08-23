#!/usr/bin/env python2
desc="""Introduce SNPs into FastA.

Heterozygous positions will be coded in small font. 
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow, 13/04/2017
"""

import os, sys, gzip
import numpy as np
from datetime import datetime
from FastaIndex import FastaIndex

#bases = "ACGTN"
bases = "ATUGCKMRYSWBVHDXN.-"
bases += bases.lower()
b2i = {b: i for i, b in enumerate(bases, 1)}
i2b = {i: b for i, b in enumerate(bases, 1)}

def fasta_streamer(handle):
    """Stream FastA lines"""
    i = 0
    name = ""
    for l in handle:
        if l.startswith(">"):
            i = 0
            name = l[1:].split()[0]
        elif l.rstrip():
            seq = l.rstrip()
            yield name, seq, i
            i += len(seq)

def load_snps(fastafn, snpfn):
    """Return chr2snps"""
    # get fasta index
    faidx = FastaIndex(fastafn)
    # generate empty arrays with 1byte per char 0-255
    chr2snps = {c: np.zeros(stats[0], dtype="uint8")
                for c, stats in faidx.id2stats.iteritems()}
    # populate with SNPs
    for l in gzip.open(snpfn):
        if l.startswith('#'): continue
        # 1       20000165        15      T       A/T     0.533/0.467
        chrom, pos, cov, ref, alt, freqs = l[:-1].split('\t')[:6]
        pos = int(pos)
        # hetero positions will be marked by lower letter
        if len(alt)>1:
            alts = alt.split('/')
            # use ref as alt in hetero position if ref is among alts
            if ref in alts:
                alt = ref.lower()
            # otherwise use most common alt
            else:
                alt = alts[0].lower()
        chr2snps[chrom][pos-1] = b2i[alt]
    return chr2snps

def get_alt_seq(seq, snps):
    """Return edited seq"""
    bases = []
    for b, bi in zip(seq, snps):
        if bi:
            if b.islower():
                b = i2b[bi].lower()
            else:
                b = i2b[bi]
        bases.append(b)
    return "".join(bases)
    
def heterozygous2fasta(handle, out, snpsfn, verbose):
    """Report to out FastA with included SNPs"""
    if verbose:
        sys.stderr.write("Loading snps...\n")
    # load snps
    chr2snps = load_snps(handle.name, snpsfn)

    if verbose:
        sys.stderr.write("Parsing FastA...\n")
    # report edited FastA
    for name, seq, s in fasta_streamer(handle):
        #if s>10000: break
        # write header
        if not s:
            sys.stderr.write(' %s      \r'%name)
            out.write(">%s\n"%name)
        # edit seqs only if any SNPs in this piece
        snps = chr2snps[name][s:s+len(seq)]
        if sum(snps):
            #pseq = seq
            seq = get_alt_seq(seq, snps)
            #if seq.upper()!=pseq.upper(): print pseq; print seq; print snps
        # write seq
        out.write(seq+"\n")
         
def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-i", "--input", required=1, type=file, help="input FastA stream")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("w"), help="output stream [stdout]")
    parser.add_argument("-s", "--snps", required=1, help="SNP files")

    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    heterozygous2fasta(o.input, o.output, o.snps, o.verbose)
        
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
