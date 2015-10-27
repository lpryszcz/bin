#!/usr/bin/env python
desc="""Report SNPs from amplicons. 
 
CHANGELOG:
v1.1
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Bratislava, 17/09/2015
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
    # count also bases like reference
    if base_ref!="N":
        basei = alphabet.index(base_ref.upper().index())
        baseCounts[basei] += alts.count('.')
        baseCounts[basei] += alts.count(',')
    #get base frequencies
    for base_count, base in sorted(baseCounts):
        freq = base_count*1.0/len(alts)#cov
        if freq >= minFreq: #base!=base_ref
            #check if alt base in both strands
            if bothStrands: 
                if not base.upper() in alts or not base.lower() in alts:
                    return
            return (base, freq) # base!=base_ref and

def get_major_alleles(base_ref, cov, alg, minFreq, alphabet, bothStrands):
    """Return major alleles that passed filtering and their frequencies."""
    #remove deletions
    alts = alg
    dels = alts.count('*') 
    #remove insertions
    alts = _remove_indels( alts )
    #get base frequencies
    bases, freqs = [], []
    # patch for ref
    alts = alts.replace(".", base_ref).replace(",", base_ref).upper()
    for base_count, base in zip((alts.count(b) for b in alphabet), alphabet):
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
    
def mpileup2calls(ref, data, minDepth, minFreq, bothStrands, \
                  alphabet, null="-"):
    """Return base calls from mpileup"""
    calls = []
    for cov, alg, qual in zip(data[0::3], data[1::3], data[2::3]):
        cov = int(cov)
        if cov<minDepth:
            calls.append(null)
            continue
        # check for SNP
        #print bases, freqs
        bases, freqs = get_major_alleles(ref, cov, alg, minFreq, alphabet, bothStrands)
        if not bases:
            calls.append(null)
            continue
        # remove ref base from alternative bases
        calls.append("".join(sorted(bases)))
    return calls
        
def parse_mpileup(bams, fasta, minDepth, minFreq, mpileup_opts, verbose, \
                  bothStrands=0, alphabet='ACGT'):
    """Run mpileup subprocess and parse output."""
    #open subprocess
    args = ['samtools', 'mpileup', "-f", fasta] + mpileup_opts.split() + bams
    if verbose:
        sys.stderr.write("Running samtools mpileup...\n %s\n" % " ".join(args))
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, bufsize=65536)
    genotyped = [0]*len(bams)
    for line in proc.stdout:
        line      = line.strip()
        lineTuple = line.split('\t')
        #get coordinate
        contig, pos, baseRef = lineTuple[:3]
        # use reference fasta
        baseRef, refFreq, refCov, gmeanQ = baseRef, set([1.0]), 100, 40.0
        baseRef = baseRef.upper()
        if baseRef not in "ACGT":
            continue
        samplesData = lineTuple[3:]
        # get calls
        calls = mpileup2calls(baseRef, samplesData, minDepth, minFreq, \
                              bothStrands, alphabet)
        uniqCalls = set(filter(lambda x: x!="-", calls))
        if len(uniqCalls)>1 or uniqCalls and baseRef not in uniqCalls:
            #print len(calls), calls
            yield tuple([contig, pos, baseRef] + calls)
        # stats
        for i, c in enumerate(calls):
            if c!="-":
                genotyped[i] += 1
    # report genotyped
    sys.stderr.write("#sample\tgenotyped positions\n")
    for b, g in zip(bams, genotyped):
        sys.stderr.write("%s\t%s\n"%(b, g))

def main():

    usage  = "%(prog)s [options] -f ref.fa -i *.bam" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.11')
    parser.add_argument("-i", "--input", nargs="+", 
                        help="input BAM file(s)")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("w"), 
                        help="output stream [stdout]")
    parser.add_argument("-f", "--fasta", required=True, 
                        help="reference FASTA file")
    parser.add_argument("--minDepth", default=5,  type=int,
                        help="minimal depth of coverage [%(default)s]")
    parser.add_argument("--minFreq",  default=0.2, type=float,
                        help="min frequency for alternative base [%(default)s]")
    parser.add_argument("--mpileup_opts",   default="-I -q 15 -Q 20",  
                        help="options passed to mpileup         [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
    
    #parse pileup in single-thread mode 
    parser = parse_mpileup(o.input, o.fasta, o.minDepth, o.minFreq, o.mpileup_opts, o.verbose)

    o.output.write("#chrom\tpos\tref\t%s\n"%"\t".join(f for f in o.input))
    info = "%s\t%s\t%s"+"\t%s"*len(o.input)+"\n"
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
