#!/usr/bin/env python
desc="""Recover telomers from short reads.
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 14/09/2016
"""

import math, os, resource, sys
from datetime import datetime
#from Bio import SeqIO
from collections import Counter
import numpy as np

nucleotides = 'ACGT'
DNAcomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

# from RapSi
def memory_usage():
    """Return memory usage in MB"""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

def encode(n, alphabet, length, base=""):
    """Converts a positive integer to a string based in given alphabet."""
    while n != 0:
        n, i = divmod(n, len(alphabet))
        base = alphabet[i] + base
    return base.rjust(length).replace(' ',alphabet[0])
        
def decode(base, alphabet, n=0):
    """Convert string based in given alphabet to int."""
    for i, char in enumerate(base, 1):
        n += alphabet.index(char) * (len(alphabet) ** (len(base) - i))
    return n  
    
def reverse_complement(mer):
    """Return DNA reverse complement"""
    return "".join(DNAcomplement[b] for b in reversed(mer))

def get_entropy(seq):
    """Return Shannon entropy of given seq"""
    M = float(len(seq))
    base2count = Counter(seq)
    entropy = -sum(c/M*(math.log(c/M, 2)) for c in base2count.itervalues())
    return entropy
    
def dnaseq2mers(seq, kmer, step, entropy=1.0, alphabet=list(nucleotides)):
    """Kmers generator for DNA seq"""
    mers = set()
    for s in xrange(0, len(seq)-kmer, step):
        mer = seq[s:s+kmer]
        # skip kmer if N inside or entropy too low        
        if "N" in mer or get_entropy(mer) < entropy:
            continue
        meri = min((decode(mer, alphabet), decode(reverse_complement(mer), alphabet)))
        mers.add(meri)
    return mers

def get_common_mers(handle, kmer, step, limit, entropy, verbose):
    """Return most common kmers"""
    mer2count = Counter()
    fq2seq    = (l.strip() for j, l in enumerate(handle) if j%4==1)
    for i, seq in enumerate(fq2seq, 1):
        if limit and i>limit:
            break
        if verbose and not i%1e3:
            sys.stderr.write(" %s %s kmers [%s MB]\r"%(i, len(mer2count), memory_usage()))
        # skip read if entropy too low
        if get_entropy(seq) < entropy:
            #print "Low entropy %s: %s" % (get_entropy(seq), seq)
            continue
        # add mers
        for m in dnaseq2mers(seq, kmer, step, entropy):
            mer2count[m] += 1
    if verbose:
        sys.stderr.write(" %s %s kmers [%s MB]\n"%(i, len(mer2count), memory_usage()))
    return mer2count
    
def fastq2telomers(handle, out, kmer, step, limit, entropy, verbose):
    """Return likely telomer sequences"""
    # get common kmers
    mer2count = get_common_mers(handle, kmer, step, limit, entropy, verbose)
    for meri, c in mer2count.most_common(50):
        mer = encode(meri, nucleotides, kmer)
        print mer, c

def main(): 
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.4a')   
    parser.add_argument("-i", "--input", default=sys.stdin, type=file, 
                        help="input stream [stdin]")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("w"), 
                        help="input stream [stdout]")
    parser.add_argument("-k", "--kmer", default=21, type=int, 
                        help="word size [%(default)s]")
    parser.add_argument("-s", "--step", default=1, type=int, 
                        help="step [%(default)s]")
    parser.add_argument("-l", "--limit", default=0, type=int, 
                        help="process only this amount of reads [all]")
    parser.add_argument("-e", "--entropy", default=1.0, type=float, 
                        help="min Shannon entropy of kmer [%(default)s]")
                        
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    # 
    fastq2telomers(o.input, o.output, o.kmer, o.step, o.limit, o.entropy, o.verbose)
    
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
