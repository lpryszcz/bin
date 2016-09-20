#!/usr/bin/env python
desc="""Find putative telomers in short reads.
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
##
    
def count_mers(handle, kmer, step, limit, entropy, verbose):
    """Return most common kmers"""
    fq2seq = (l.strip() for j, l in enumerate(handle) if j%4==1)
    mer2count = Counter()
    for i, seq in enumerate(fq2seq, 1):
        # consider adding mer2count limit instead of read limit #len(mer2count)>1e5 or 
        if limit and len(mer2count) > limit:
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

def fw(mer, alphabet='ACGT'):
    """Generator of forward mers for given mer"""
    for base in alphabet:
        yield mer[1:] + base

def bw(mer, alphabet='ACGT'):
    """Generator of reverse mers for given mer"""
    for base in alphabet:
        yield base + mer[:-1]
    
def contig2string(mers):
    """Return contig sequence, this is a concatenation of first mer and
    last character of following mers
    """
    return mers[0] + ''.join(mer[-1] for mer in mers[1:])

def get_contig_forward(mer2count, mer):
    """Return mers of contig that given mer belongs to"""
    mers = [mer]
    while True:
        # if forward-degree of the last mer (sum of True/False) != 1 then break the contig
        if sum(x in mer2count for x in fw(mers[-1])) != 1:
            break
        _mer = [x for x in fw(mers[-1]) if x in mer2count][0]
        # break out of cycles or mobius contigs OR harpins OR not forward mer
        if _mer == mer or _mer == twin(mer) or _mer == twin(mers[-1]) or \
           sum(x in mer2count for x in bw(_mer)) != 1:
            break
        # store next mers
        mers.append(_mer)
    return mers
    
def get_contig(mer2count, mer):
    """Return sequence of contig and list of its kmers"""
    mers = get_contig_forward(mer2count, mer)
    # get backward contigs if not forward
    if mer not in fw(mers[-1]):
        c_bw = get_contig_forward(mer2count, twin(mer))
        mers = [twin(x) for x in c_bw[-1:0:-1]] + mers
    return contig2string(mers), mers
    
def get_contigs(mer2count):
    """Return sequences of contigs
    Adapted from: https://github.com/pmelsted/dbg
    """
    contigs, coverages, done = [], [], set()
    for mer in mer2count:
        # skip already processed
        if mer in done:
            continue
        # get contig and its kmers
        seq, mers = get_contig(mer2count, mer)
        contigs.append(seq)        
        # get kmer coverage
        cov = 1.0 * sum(mer2count[mer] for mer in mers) / len(mers)
        coverages.append(cov)
        # update done
        done.update(mers + map(reverse_complement, mers))
    return contigs, coverages

def getsubs(loc, s):
    """Get substring"""
    substr = s[loc:]
    i = -1
    while(substr):
        yield substr
        substr = s[loc:i]
        i -= 1

def get_longest_repetitive_substring(seq, minocc=2):
    """Return longest repetitive substring and its number of occurrences.
    Return empty str and 0 if no substring having at least minocc occurrences is found. 
    Adapted from http://stackoverflow.com/a/11091454/632242
    """
    occ = defaultdict(int)
    # tally all occurrences of all substrings
    for i in range(len(seq)):
        for sub in getsubs(i, seq):
            occ[sub] += 1

    # filter out all substrings with fewer than minocc occurrences
    occ_minocc = [k for k, v in occ.items() if v >= minocc]
    # return longest 
    if occ_minocc:
        maxkey =  max(occ_minocc, key=len)
        return maxkey, occ[maxkey]
    return '', 0
    
def get_telomers(seq, minlength=5):
    """Return lists of the longest repeatitive substrings and likely telomeric repeat"""
    repeats = []
    # get longest repetitive substrings
    _seq = seq
    while _seq: 
        _seq, i = get_longest_repetitive_substring(_seq)
        if len(_seq) > minlength:
            repeats.append(_seq)
    # check telomers
    if repeats:
        repeat = repeats[-1]
        starts = [m.start() for m in re.compile(repeat).finditer(seq)]
        if len(starts)<2:
            return repeats, ''
        trs = []
        for i in range(len(starts)-1):
            s, e = starts[i:i+2]
            tr = seq[s:e]
            # check if continues
            s1, s2 = tr, seq[e:e+len(tr)]
            ## check if telomeric repeat is before and after - this is quite restrictive!
            # make sure s1 is shorter than s2 (contig may be short)
            if len(s2)<len(s1):
                s1, s2 = s2, s1
            e1, e2 = tr, seq[0:s]
            if len(e2)<len(e1):
                e1, e2 = e2, e1
            # only valid telomeric repeat if telomeric repeat sequence is repeated
            if s2.startswith(s1) and e2.endswith(e1):
                trs.append(tr)
            else:
                return repeats, ''
        if trs and len(set(trs))==1:
            tr = trs[0]
            return repeats, tr
    return repeats, ''
    
def fastq2telomers(handle, out, kmer, step, limit, minlength, topmers, entropy, verbose):
    """Return likely telomer sequences"""
    # get common kmers
    if verbose:
        sys.stderr.write("Counting kmers in reads...\n")
    mer2count, reads = count_mers(handle, kmer, step, limit, entropy, verbose)
    if reads < 1e5:
        sys.stderr.write("[WARNING] Only %s reads processed! Consider processing large pool of reads...\n"%reads)

    meri = sum(mer2count.itervalues())
    if verbose: 
        sys.stderr.write("%s occurrences of %s %s-kmers\n"%(meri, len(mer2count), kmer))
        #for mer, c in mer2count.most_common(10):
        #    sys.stderr.write(" %s %s\n" % (mer, c))

    # remove rare kmers
    mincount = mer2count.most_common(topmers)[-1][1]
    mincount = mincount if mincount > 1 else 2
    if verbose:
        sys.stderr.write("Removing rare kmers (< %s occurrences) ...\n"%mincount)
    mer2count = {mer: c for mer, c in mer2count.iteritems() if c >= mincount}
    if verbose:
        sys.stderr.write(" %s kmers kept.\n"%len(mer2count))
        
    # assemble common kmers
    if verbose:
        sys.stderr.write("Assembly...\n")
    contigs, coverages = get_contigs(mer2count)
    if verbose:
        sys.stderr.write(" %s bp in %s contigs [%s Mb].\n"%(sum(len(c) for c in contigs), len(contigs), memory_usage()))
        
    # select most likely telomere
    if verbose:
        sys.stderr.write("Reporting most likely telomeric repeats...\n")
    k = 0
    subrepeats, trs = [], []
    for seq in contigs:
        repeats, tr = get_telomers(seq, minlength)
        subrepeats.append(repeats)
        trs.append(tr)
        if tr:
            k +=1
    i = 0
    for i, (cov, seq, repeats, tr) in enumerate(sorted(zip(coverages, contigs, subrepeats, trs), reverse=1), 1):
        if i==1:
            maxcov = cov
        # report only trs if any present
        if k and not tr:
            continue
        repeats, telomeric_repeat = get_telomers(seq, minlength)
        header = "contig%d cov:%.2f telomeric_repeat:%s|%s substrings:%s"%(i, cov, telomeric_repeat, reverse_complement(telomeric_repeat), ";".join(repeats), )
        out.write('>%s\n%s\n'%(header, seq))
    if verbose and i:
        sys.stderr.write(" %s contigs & %s putative telomers stored.\n"%(i, k))
        sys.stderr.write("  contigs kmer coverage: %.2f - %.2f.\n"%(cov, maxcov))
        
def main(): 
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0a')   
    parser.add_argument("-i", "--input", default=sys.stdin, type=file, 
                        help="input stream [stdin]")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("w"), 
                        help="input stream [stdout]")
    parser.add_argument("-k", "--kmer", default=41, type=int, 
                        help="word size [%(default)s]")
    parser.add_argument("-s", "--step", default=1, type=int, 
                        help="step [%(default)s]")
    parser.add_argument("-l", "--limit", default=1e6, type=float, 
                        help="process until reaching this amount of kmers [%(default)s]")
    parser.add_argument("-e", "--entropy", default=0.75, type=float, 
                        help="min Shannon entropy of kmer [%(default)s]")
    #parser.add_argument("-f", "--minfreq", default=1e-05, type=float, 
    #                    help="min frequency of kmer [%(default)s]")
    parser.add_argument("-m", "--minlength", default=5, type=int, 
                        help="min telomer length [%(default)s]")
    parser.add_argument("-t", "--topmers", default=1000, type=int, 
                        help="no. of top occurring kmer to use [%(default)s]")
                        
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    # 
    fastq2telomers(o.input, o.output, o.kmer, o.step, o.limit, o.minlength,
                   o.topmers, o.entropy, o.verbose)
    
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
