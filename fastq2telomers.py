#!/usr/bin/env python
desc="""Recover telomers from short reads.

- mer as int 1:54 104Mb
- mer as str 0:56 140Mb
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
twin = reverse_complement

def get_entropy(seq):
    """Return Shannon entropy of given seq"""
    M = float(len(seq))
    base2count = Counter(seq)
    entropy = -sum(c/M*(math.log(c/M, 2)) for c in base2count.itervalues())
    return entropy
    
def dnaseq2mers(seq, kmer, step, entropy=1.0, alphabet=list(nucleotides)):
    """Kmers generator for DNA seq"""
    mers = [] #set()
    for _seq in filter(lambda x: len(x)>=kmer, seq.split('N')):
        for s in xrange(0, len(_seq)-kmer, step):
            mer = seq[s:s+kmer]
            # skip kmer if N inside or entropy too low        
            if get_entropy(mer) < entropy:
                continue
            # store kmer and its reverse complement
            mers += [mer, reverse_complement(mer)]
    return mers

def count_mers(handle, kmer, step, limit, entropy, verbose):
    """Return most common kmers"""
    fq2seq = (l.strip() for j, l in enumerate(handle) if j%4==1)
    mer2count = Counter()
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

## assembly
def fw(km):
    for x in 'ACGT':
        yield km[1:]+x

def bw(km):
    for x in 'ACGT':
        yield x + km[:-1]
    
def contig_to_string(c):
    return c[0] + ''.join(x[-1] for x in c[1:])

def get_contig(d,km):
    c_fw = get_contig_forward(d,km)
    
    c_bw = get_contig_forward(d,twin(km))

    if km in fw(c_fw[-1]):
        c = c_fw
    else:
        c = [twin(x) for x in c_bw[-1:0:-1]] + c_fw
    return contig_to_string(c),c
        

def get_contig_forward(d,km):
    c_fw = [km]
    
    while True:
        if sum(x in d for x in fw(c_fw[-1])) != 1:
            break
        
        cand = [x for x in fw(c_fw[-1]) if x in d][0]
        if cand == km or cand == twin(km):
            break # break out of cycles or mobius contigs
        if cand == twin(c_fw[-1]):
            break # break out of hairpins
        
        if sum(x in d for x in bw(cand)) != 1:
            break

        c_fw.append(cand)

    return c_fw

def all_contigs(d,k):
    done = set()
    r = []
    for x in d:
        if x not in done:
            s,c = get_contig(d,x)
            for y in c:
                done.add(y)
                done.add(twin(y))
            r.append(s)
    
    G = {}
    heads = {}
    tails = {}
    for i,x in enumerate(r):
        G[i] = ([],[])
        heads[x[:k]] = (i,'+')
        tails[twin(x[-k:])] = (i,'-')
    
    for i in G:
        x = r[i]
        for y in fw(x[-k:]):
            if y in heads:
                G[i][0].append(heads[y])
            if y in tails:
                G[i][0].append(tails[y])
        for z in fw(twin(x[:k])):
            if z in heads:
                G[i][1].append(heads[z])
            if z in tails:
                G[i][1].append(tails[z])

    return G,r

def print_GFA(G,cs,k):
    print "H  VN:Z:1.0"
    for i,x in enumerate(cs):
        print "S\t%d\t%s\t*"%(i,x)
        
    for i in G:
        for j,o in G[i][0]:
            print "L\t%d\t+\t%d\t%s\t%dM"%(i,j,o,k-1)
        for j,o in G[i][1]:
            print "L\t%d\t-\t%d\t%s\t%dM"%(i,j,o,k-1)

def save_contigs(out, cs):
    for i, x in enumerate(cs, 1):
        out.write('>contig%d\n%s\n'%(i, x))
##
    
def fastq2telomers(handle, out, kmer, step, limit, entropy, verbose):
    """Return likely telomer sequences"""
    # get common kmers
    mer2count = count_mers(handle, kmer, step, limit, entropy, verbose)
    for mer, c in mer2count.most_common(10):
        print mer, c

    # remove not frequent kmers
    for mer in filter(lambda x: mer2count[x]<10, mer2count):
        del mer2count[mer]
        
    # assemble common kmers
    G, cs = all_contigs(mer2count, kmer)
    save_contigs(out, cs) 

    # select most likely telomere - the contig having repeats
    

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
    parser.add_argument("-l", "--limit", default=0, type=float, 
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
