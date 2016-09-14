#!/usr/bin/env python
desc="""Recover telomers from short reads.

- mer as int 1:54 104Mb
- mer as str 0:56 140Mb
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 14/09/2016
"""

import math, os, re, resource, sys
from datetime import datetime
from collections import Counter, defaultdict

nucleotides = 'ACGT'
DNAcomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

npat = re.compile(r'N+')

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
    mers = [] 
    for _seq in filter(lambda x: len(x)>=kmer, npat.split(seq)):
        for s in range(0, len(_seq)-kmer, step):
            mer = _seq[s:s+kmer]
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
            continue
        # add mers
        for m in dnaseq2mers(seq, kmer, step, entropy):
            mer2count[m] += 1
    if verbose:
        sys.stderr.write(" %s %s kmers [%s MB]\n"%(i, len(mer2count), memory_usage()))
    return mer2count, i

## de Bruijn assembly
# adapted from https://pmelsted.wordpress.com/2013/11/23/naive-python-implementation-of-a-de-bruijn-graph/
def fw(km):
    for x in 'ACGT':
        yield km[1:]+x

def bw(km):
    for x in 'ACGT':
        yield x + km[:-1]
    
def contig2string(c):
    return c[0] + ''.join(x[-1] for x in c[1:])

def get_contig_forward(mer2count, mer):
    c_fw = [mer]
    while True:
        if sum(x in mer2count for x in fw(c_fw[-1])) != 1:
            break
        cand = [x for x in fw(c_fw[-1]) if x in mer2count][0]
        # break out of cycles or mobius contigs
        if cand == mer or cand == twin(mer):
            break
        # break out of hairpins
        if cand == twin(c_fw[-1]):
            break 
        if sum(x in mer2count for x in bw(cand)) != 1:
            break
        c_fw.append(cand)
    return c_fw

def get_contig(mer2count, mer):
    """Return sequence of contig and list of its kmers"""
    c_fw = get_contig_forward(mer2count, mer)
    if mer in fw(c_fw[-1]):
        c = c_fw
    else:
        c_bw = get_contig_forward(mer2count, twin(mer))
        c = [twin(x) for x in c_bw[-1:0:-1]] + c_fw
    return contig2string(c), c

def get_contigs(mer2count):
    """Return sequences of contigs"""
    contigs, done = [], set()
    for mer in mer2count:
        if mer in done:
            continue
        # get contig and its kmers
        seq, mers = get_contig(mer2count, mer)
        contigs.append(seq)
        # update done
        done.update(mers + map(reverse_complement, mers))
    return contigs

##
# adapted from http://stackoverflow.com/a/11091454/632242
def getsubs(loc, s):
    substr = s[loc:]
    i = -1
    while(substr):
        yield substr
        substr = s[loc:i]
        i -= 1

def longestRepetitiveSubstring(r, minocc=2):
    occ = defaultdict(int)
    # tally all occurrences of all substrings
    for i in range(len(r)):
        for sub in getsubs(i,r):
            occ[sub] += 1

    # filter out all substrings with fewer than minocc occurrences
    occ_minocc = [k for k,v in occ.items() if v >= minocc]

    if occ_minocc:
        maxkey =  max(occ_minocc, key=len)
        return maxkey, occ[maxkey]
    #else:
    #    raise ValueError("no repetitions of any substring of '%s' with %d or more occurrences" % (r, minocc))
    return '', 0
##
    
def get_telomers(seq, minlength=5):
    """Return likely telomers or telomer repeats"""
    telomers = []
    while seq: 
        seq, i = longestRepetitiveSubstring(seq)
        if len(seq) > minlength:
            telomers.append(seq)
    return telomers
    
def fastq2telomers(handle, out, kmer, step, limit, minlength, minfreq, entropy, verbose):
    """Return likely telomer sequences"""
    # get common kmers
    if verbose:
        sys.stderr.write("Counting kmers in reads...\n")
    mer2count, reads = count_mers(handle, kmer, step, limit, entropy, verbose)
    if reads < 1e5:
        sys.stderr.write("[WARNING] Only %s reads processed! Consider processing large pool of reads...\n"%reads)

    meri = sum(mer2count.itervalues())
    if verbose:
        sys.stderr.write("%s occurences of %s %s-kmers\n"%(meri, len(mer2count), kmer))
        for mer, c in mer2count.most_common(10):
            sys.stderr.write(" %s %s\n" % (mer, c))

    # remove rare kmers
    mincount = round(minfreq*meri if minfreq*meri > 1 else 1)
    if verbose:
        sys.stderr.write("Removing rare kmers (< %s) ...\n"%minfreq)
    for mer in filter(lambda x: mer2count[x] <= mincount, mer2count):
        del mer2count[mer]
    if verbose:
        sys.stderr.write(" %s kmers kept.\n"%len(mer2count))
        
    # assemble common kmers
    if verbose:
        sys.stderr.write("Assembly...\n")
    contigs = get_contigs(mer2count)
        
    # select most likely telomere
    if verbose:
        sys.stderr.write("Reporting most likely telomeres...\n")
    k = 0
    for i, seq in enumerate(contigs, 1):
        telomers = get_telomers(seq, minlength)
        k += len(telomers)
        out.write('>contig%d putative telomers: %s\n%s\n'%(i, "; ".join(telomers), seq))
    if verbose:
        sys.stderr.write(" %s contigs & %s putative telomers stored.\n"%(i, k))
        
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
    parser.add_argument("-k", "--kmer", default=41, type=int, 
                        help="word size [%(default)s]")
    parser.add_argument("-s", "--step", default=1, type=int, 
                        help="step [%(default)s]")
    parser.add_argument("-l", "--limit", default=0, type=float, 
                        help="process only this amount of reads [all]")
    parser.add_argument("-e", "--entropy", default=1.0, type=float, 
                        help="min Shannon entropy of kmer [%(default)s]")
    parser.add_argument("-f", "--minfreq", default=1e-05, type=float, 
                        help="min frequency of kmer [%(default)s]")
    parser.add_argument("-m", "--minlength", default=5, type=int, 
                        help="min telomer length [%(default)s]")
                        
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    # 
    fastq2telomers(o.input, o.output, o.kmer, o.step, o.limit, o.minlength,
                   o.minfreq, o.entropy, o.verbose)
    
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
