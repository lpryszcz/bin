#!/usr/bin/env python
desc="""Count k-mers occurencies in reads and report freq plot.
Not working for k > 16...

NOT RECOMMENDED! USE JELLYFISH INSTEAD:
#jellyfish
k=31; q=20; 
for f in _archives_CANTHE/CAN*read?.fastq.gz; do 
 s=`echo $f | cut -f2 -d"/" | cut -f1 -d"."`; 
 if [ ! -s _archives/$s.k$k.q$q.hist ]; then 
  echo `date` $f $s; 
  zcat $f | fastq2fasta.py -v -l $k -q $q | jellyfish count -m31 -s1000000000 -o $s.k$k.q$q.jelly /dev/stdin; jellyfish histo -h10000 _archives/$s.k$k.q$q.jelly_0 > _archives/$s.k$k.q$q.hist
 fi; done
plot_2d.py -v -i _archives/*.k31.q20.hist -o CANME.k31.q20.png


ver 1.1:
- quality control
- numpy.array for counting
- multiprocessing
ver 1.0:
- Biopython FASTQ parser
- 
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona/Mizerow, 9/05/2013
"""

import argparse, gzip, os, resource, sys
from datetime import datetime
from Bio      import SeqIO
from matplotlib import pyplot as plt
from multiprocessing import Pool
import numpy as np

def get_mer(mer):
    """Return forward or reverse sequence"""
    merC = mer.complement()
    for a, b in zip(mer, merC):
        if   a<b:
            return mer
        elif b<a:
            return merC
    return mer

def init_args(*args):
    global k, ns, base2digit
    k, ns, base2digit = args
    
def seq2mers(r):
    """Return list of kmer ids from given read sequence"""
    ids = []
    for j in range(len(r.seq)-k):
        mer = get_mer(r.seq[j:j+k])
        if not ns and "N" in mer:
            continue
        #get id and store it
        mid = int("".join(base2digit[base] for base in mer), len(base2digit))
        #if mid < len(base2digit)**k/2:
        ids.append(mid)
    return ids
    
def get_kmer_counts(input, output, k, ns, nprocs, verbose):
    """Analyse kmers. Multiprocessing enabled"""
    #define base2digit dict for 4-char seq
    base2digit = {"A": "0", "C": "1", "G": "2", "T": "3"}    
    if ns:
        #change to 5-char seq if Ns in seq
        base2digit = {"A": "0", "C": "1", "G": "2", "N": "3", "T": "4"}
    #init mer counts
    #255 for uint8 #65,535 for uint16 or #4,294,967,295 for uint32 
    merCounts = np.zeros(len(base2digit)**k/2, dtype='uint16')
    #start pool #maxtasksperchild=1000)
    p = Pool(nprocs, initializer=init_args, initargs=(k, ns, base2digit)) 
    #process reads
    for i, ids in enumerate(p.imap_unordered(seq2mers, SeqIO.parse(input, 'fastq'), \
                                             chunksize=100), 1):
        if not i%1e4:
            sys.stderr.write(" %s [%s Mb]\r"%(i, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
        for mid in ids:
            merCounts[mid] += 1
    sys.stderr.write(" %s [%s Mb]\n"%(i, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
    #get mer freq
    maxCount    = merCounts.max()
    if maxCount < 100:
        maxCount = 100
    occurencies = [0]*maxCount
    for c in merCounts:
        occurencies[c-1] += 1
    #write to file
    output.write("\n".join("%s\t%s"%xy for xy in enumerate(occurencies,1))+"\n")
    return occurencies
    
'''v1.0    
def get_kmer2count(input, kmer, ns, verbose):
    """Count kmer occurencies. Uses 2 threads:
    - fastq parsing
    - kmers storing
    """
    kmer2count = {}
    step = 10**3
    pI   = 0
    #parse fastq
    for i, r in enumerate(SeqIO.parse(input,'fastq'),1):
        #process mers
        for j in range(len(r.seq)-kmer):
            merF = str(r.seq[j:j+kmer])
            if not ns and "N" in merF:
                continue
            merR = str(r.seq[j:j+kmer].reverse_complement())        
            #add to dict
            if   merF in kmer2count:
                kmer2count[merF] += 1
            elif merR in kmer2count:
                kmer2count[merR] += 1
            else:
                kmer2count[merF]  = 1 
        #print info        
        if i>pI:
            pI += step
            sys.stderr.write(" %12i [%s Mb]\r" % (i, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
    sys.stderr.write(" %12i [%s Mb] kmer2count: %s \n" % (i, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, len(kmer2count)))
    return kmer2count

def get_occurencies(kmer2count, output, verbose):
    """ """
    #max kmer occurencies
    maxCount    = max(c for k,c in kmer2count.items())
    if maxCount < 100:
        maxCount = 100
    occurencies = [0]*maxCount
    for k,c in kmer2count.items():
        occurencies[c-1] += 1
    #write to file
    output.write("\n".join("%s\t%s"%xy for xy in enumerate(occurencies,1))+"\n")
    return occurencies

def seq2int(seq, minseqi, base2digit):
    """Return int representation of seq."""
    #convert seq into seq of digits
    dseq = "".join(base2digit[base] for base in seq)
    #if not Ns in seq, encode as 4-char int
    seqi = int("1"+dseq, len(base2digit))
    #return iseq minus int for min seq, that is A*kmer
    return seqi-minseqi
        
def get_kmer_counts(input, kmer, ns, verbose):
    """Count kmer occurencies.
    sys.getsizeof([0]*4**15) is just 1,073,741,824 elements 10e9
    while 4**31=10e18!
    Out[35]: 8,589,934,664 -> 8.5Gb for k=15
    use numpy
    """
    #define base2digit dict for 4-char seq
    base2digit = {"A": "0", "C": "1", "G": "2", "T": "3"}    
    if ns:
        #change to 5-char seq if Ns in seq
        base2digit = {"A": "0", "C": "1", "G": "2", "N": "3", "T": "4"}
    #define lowest & highest possible int representation of seq (AAAAA)
    ##easier would be int("1"+"0"*kmer,4) and int("3"*kmer,4)
    minseqi = seq2int("A"*kmer,0,base2digit) 
    maxseqi = seq2int("T"*kmer,minseqi,base2digit)
    #print minseqi, maxseqi
    #define empty occurencies list
    counts = np.zeros((maxseqi+1)/2, dtype=int) #[0] * maxseqi
    #parse fastq
    step   = 10**3
    pI     = 0
    for i,r in enumerate(SeqIO.parse(input,'fastq'),1):
        #process mers
        for j in range(len(r.seq)-kmer):
            mer  = str(r.seq[j:j+kmer]).upper()
            if not ns and "N" in mer:
                continue
            merR = str(r.seq[j:j+kmer].reverse_complement())    
            #get smaller seq
            if   merR < mer:
                mer = merR
            #update count
            meri = seq2int(mer, minseqi, base2digit)
            counts[meri] += 1
        #print info        
        if i>pI:
            pI += step
            sys.stderr.write(" %12i [%s Mb]\r" % (i, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
    sys.stderr.write(" %12i [%s Mb] kmers: %s \n" % (i, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, len(kmer2count)))
    return counts
'''
def plot(occurencies, kmer, output, title):
    """Plot figure."""
    #get figure
    plt.figure(figsize=(6,6))
    plt.plot(xrange(1,101),occurencies[:100],linewidth=2.0)
    #add title and axis labels
    plt.title(title+" [k=%s]"%kmer)
    plt.xlabel('k-mer frequency')
    plt.ylabel('number of k-mers with this frequency')        
    #show plot if not outfile provided
    if output.name=="<stdout>":
        plt.show()
    else:
        format = "png"
        fpath  = "%s.%s" % (output.name, format)
        plt.savefig(fpath, dpi=300, facecolor='w', edgecolor='w',\
          orientation='landscape', format=format, transparent=False)
        
def fastq2kmers(input, output, kmer, ns, nprocs, title, verbose):        
    """Process fastq from stdin."""
    if verbose:
        sys.stdout.write("[fastq2kmers] Parsing fastq...\n")
    occurencies = get_kmer_counts(input, output, kmer, ns, nprocs, verbose)
    '''
    kmer2count  = get_kmer2count(input, kmer, ns, verbose)
    #get occuriences
    if verbose:
        sys.stdout.write("[fastq2kmers] Getting frequencies...\n")
    occurencies = get_occurencies(kmer2count, output, verbose)'''

    #plot
    if verbose:
        sys.stdout.write("[fastq2kmers] Plotting...\n")
    plot(occurencies, kmer, output, title)
    
def main():
    usage   = "zcat fastq.gz | %(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",   default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-i", dest="input",     default=sys.stdin, type=file, 
                        help="input stream      [stdin]")
    parser.add_argument("-o", dest="output",    default=sys.stdout, type=argparse.FileType("w"),
                        help="output stream     [stdout]")
    parser.add_argument("-k", dest="kmer",      default=31, type=int, 
                        help="k-mer             [%(default)s]")
    parser.add_argument("-t", dest="title",     default="",
                        help="plot title        [%(default)s]")
    parser.add_argument("-n", "--nprocs",         default=4,
                        help="no. of processes  [%(default)s]")
    parser.add_argument("--include-Ns", dest="ns", default=False, action="store_true", 
                        help="include k-mers with Ns [%(default)s]")
                         
    o = parser.parse_args()
    if o.verbose:
        sys.stdout.write( "[fastq2kmers] Options: %s\n" % str(o) )

    fastq2kmers(o.input, o.output, o.kmer, o.ns, o.nprocs, o.title, o.verbose)
	
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\n[fastq2kmers] Ctrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stdout.write( "[fastq2kmers] Time elapsed: %s\n" % dt )
