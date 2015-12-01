#!/usr/bin/env python
desc="""Identify RNA editing sites. 
 
CHANGELOG:
v1.12
- skip alt base calling if less than 3 reads (--minAltReads)
- -o / --output option added
v1.11
- use fasta file as reference

TBD:
- rescue major alt haplotype if the other is only very low freq
- multithreading support
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw, 21/07/2015
"""

import argparse, os, re, sys
from datetime import datetime
import subprocess
import numpy as np

#find stats of the reads in mpileup
##http://samtools.sourceforge.net/pileup.shtml
read_start_pat = re.compile('\^.')
indel_pat = re.compile('[+-]\d+')

def _remove_indels(alts):
    """
    Remove indels from mpileup.     .$....,,,,....,.,,..,,.,.,,,,,,,....,.,...,.,.,....,,,........,.A.,...,,......^0.^+.^$.^0.^8.^F.^].^],
    ........,.-25ATCTGGTGGTTGGGATGTTGCCGCT..
    """
    #But first strip start/end read info.
    alts = "".join(read_start_pat.split(alts)).replace('$', '')
    #remove indels info
    m = indel_pat.search(alts)
    while m:
        #remove indel
        pos = m.end() + int(m.group()[1:])
        alts = alts[:m.start()] + alts[pos:]
        #get next match
        m = indel_pat.search(alts, m.start())
    return alts

def get_alt_allele(base_ref, cov, alg, minFreq, alphabet, reference, bothStrands):
    """Return alternative allele only if different than ref and freq >= minFreq."""
    #remove deletions
    alts = alg
    dels = alts.count('*') 
    #remove insertions
    alts = _remove_indels(alts)
    #get base counts
    baseCounts = [(alts.upper().count(base), base) for base in alphabet]
    #get base frequencies
    for base_count, base in sorted(baseCounts):
        freq = base_count*1.0/len(alts)#cov
        if base!=base_ref and freq >= minFreq:
            #check if alt base in both strands
            if bothStrands: 
                if not base.upper() in alts or not base.lower() in alts:
                    return
            return (base, freq) # base!=base_ref and

def get_major_alleles(cov, alg, minFreq, alphabet, bothStrands, minCount=3):
    """Return major alleles that passed filtering and their frequencies."""
    #remove deletions
    alts = alg
    dels = alts.count('*') 
    #remove insertions
    alts = _remove_indels( alts )
    #get base frequencies
    bases, freqs = [], []
    for base_count, base in zip((alts.upper().count(b) for b in alphabet), alphabet):
        # skip if alt base calling if less than 3 reads
        if base_count < minCount:
            continue
        freq = base_count*1.0/cov
        #check freq
        if freq < minFreq:
            continue
        #check if alt base in both strands
        if bothStrands: 
            if not base.upper() in alts or not base.lower() in alts:
                continue
        bases.append(base)
        freqs.append(freq)
    return bases, freqs

def concatenate_mpileup(data):
    """Concatenate mpileup from multiple samples / replicas"""
    cov, algs, quals = 0, "", ""
    for _cov, _algs, _quals in zip(data[0::3], data[1::3], data[2::3]):
        cov += int(_cov)
        if int(_cov):
            algs  += _algs
            quals += _quals
    return cov, algs, quals

def get_meanQ(quals, offset=33):
    """Return mean quality"""
    return np.mean([ord(q)-offset for q in quals])
    
def parse_mpileup(dna, rna, fasta, minDepth, minDNAfreq, minRNAfreq,
                  mpileup_opts, verbose, bothStrands=0, alphabet='ACGT'):
    """Run mpileup subprocess and parse output."""
    #open subprocess
    args = ['samtools', 'mpileup']
    if not dna:
        args += ["-f", fasta]
    args += mpileup_opts.split() + dna + rna 
    if verbose:
        sys.stderr.write("Running samtools mpileup...\n %s\n" % " ".join(args))
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, bufsize=65536)
    for line in proc.stdout:
        line      = line.strip()
        lineTuple = line.split('\t')
        #get coordinate
        contig, pos, baseRef = lineTuple[:3]
        # DNA bam files
        if dna:
            refData = lineTuple[3:3+3*len(dna)]
            samplesData = lineTuple[3+3*len(dna):]

            refCov, refAlgs, refQuals = concatenate_mpileup(refData)
            if refCov < minDepth:
                continue
            baseRef, refFreq = get_major_alleles(refCov, refAlgs, minDNAfreq, alphabet, bothStrands)
            if not baseRef:
                continue
            # mean qual
            gmeanQ = get_meanQ(refQuals)
        # use reference fasta
        else:
            baseRef, refFreq, refCov, gmeanQ = baseRef, set([1.0]), 100, 40.0
            baseRef = baseRef.upper()
            if baseRef not in "ACGT":
                continue
            samplesData = lineTuple[3:]
            
        cov, alg, quals = concatenate_mpileup(samplesData) #in zip(outFiles, samplesData[0::3], samplesData[1::3], samplesData[2::3]):
        cov = int(cov)
        if cov<minDepth:
            continue
        # check for SNP
        bases, freqs = get_major_alleles(cov, alg, minRNAfreq, alphabet, bothStrands)
        if not bases or bases==baseRef:
            continue
        # remove ref base from alternative bases
        refBase = baseRef[0]
        if refBase in bases:
            idx = bases.index(refBase)
            freqs = freqs[:idx] + freqs[idx+1:]
            bases.remove(refBase)
        #print refCov, refAlgs, refQuals, refData

        # rescue major alt haplotype if the other is only very low freq
        
        #print cov, alg, quals, samplesData
        if len(baseRef) != 1 or len(bases) != 1:
            info = "[WARNING] Wrong number of bases: %s:%s %s %s\n"
            sys.stderr.write(info%(contig, pos, ",".join(baseRef), ",".join(bases)))
            continue
        # skip if the same base
        #if not bases.difference(baseRef):
        #    continue
        meanQ = get_meanQ(quals)
        yield (contig, pos, refBase, bases.pop(), refCov, gmeanQ, refFreq.pop(), \
               cov, meanQ, freqs.pop())

def init_args(*args):
    global c2i, bufferSize, ivals, maxp, pref
    c2i, bufferSize = args
    # keep info about intervals, max position and current reference
    ivals, maxp, pref = [], 0, 0
    
def worker(args):
    """Count overlapping intervals with given read alignment.
    The algorithm support spliced alignments. """
    global c2i, bufferSize, ivals, maxp, pref
        
def main():

    usage  = "%(prog)s [options] -b ref.bam *.bam" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.12')
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-r", "--rna", nargs="+", 
                        help="input RNA-Seq BAM file(s)")
    refpar = parser.add_mutually_exclusive_group(required=True)
    refpar.add_argument("-d", "--dna", nargs="*", default = [], 
                        help="input DNA-Seq BAM file(s)")
    refpar.add_argument("-f", "--fasta", default = None, 
                        help="reference FASTA file")
    parser.add_argument("--minDepth", default=5,  type=int,
                        help="minimal depth of coverage [%(default)s]")
    parser.add_argument("--minAltReads", default=3,  type=int,
                        help="minimum no. of reads with alternative base to call RNA editing [%(default)s]")
    parser.add_argument("--minRNAfreq",  default=0.01, type=float,
                        help="min frequency for RNA editing base [%(default)s]")
    parser.add_argument("--minDNAfreq",  default=0.99, type=float,
                        help="min frequency for genomic base [%(default)s]")
    parser.add_argument("--mpileup_opts",   default="-I -q 15 -Q 20",  
                        help="options passed to mpileup         [%(default)s]")
    parser.add_argument("-t", "--threads", default=1,  type=int,
                        help="number of cores to use [%(default)s] NOT IMPLEMENTED YET!")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    #parse pileup in single-thread mode 
    if o.threads<2:
        parser = parse_mpileup(o.dna, o.rna, o.fasta, o.minDepth, o.minDNAfreq, o.minRNAfreq, \
                               o.mpileup_opts, o.verbose)
    else:
        # get chromosomes
        
        # process by chromosome
        p = Pool(o.threads, initializer=init_args, initargs=(c2i, bufferSize))
        parser = Pool.imap(worker, iterator, chunksize=10000)

    info = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
    for data in parser:
        o.output.write(info%data)
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
