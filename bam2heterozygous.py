#!/usr/bin/env python
desc="""Identify RNA editing sites from RNAseq and DNAseq alignements (.bam).
Alternatively, reference genome can be used instead of DNAseq,
but at the cost of higher false positive. 

TBD:
- editing from heterozygous sites?
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw/Bratislava/Fribourg, 21/07/2015
"""

import os, sys, pysam, resource
from datetime import datetime
from multiprocessing import Pool
import numpy as np

alphabet = "ACGT"
base2index = {b: i for i, b in enumerate(alphabet)}
for i, b in enumerate(alphabet.lower()):
    base2index[b] = i

def is_duplicate(a, pa):
    """Return True if read is duplicate"""
    if pa and a.pos==pa.pos and a.flag==pa.flag and a.isize==pa.isize \
       and a.cigarstring==pa.cigarstring and a.seq==pa.seq:
        return True
    
def is_qcfail(a, mapq=15):
    """Return True if alignment record fails quality checks"""
    if a.mapq<mapq or a.flag&3840: # a.is_duplicate or a.is_secondary or a.is_qcfail or a.is_supplementary:
        return True
        
def _match(refi, readi, bases): return refi+bases, readi+bases, True
def _insertion(refi, readi, bases): return refi, readi+bases, []
def _deletion(refi, readi, bases): return refi+bases, readi, []
code2function = {0: _match, 7: _match, 8: _match, 1: _insertion, 6: _insertion,
                 2: _deletion, 3: _deletion, 4: _insertion, 5: _insertion}

def get_blocks(a, start, end, baseq, i, basesize):
    """Return tuple of aligned position of query and reference"""
    readi, refi = 0, a.pos
    for code, bases in a.cigar:
        prefi, preadi = refi, readi
        refi, readi, data = code2function[code](refi, readi, bases)
        if not data or refi<start-1:
            continue
        if prefi<start:
            bases -= start-prefi
            preadi += start-prefi
            prefi = start
        if refi>end:
            bases -= refi-end
        if bases<1:
            break
        block = [0]*basesize*bases 
        for ii, (b, q) in enumerate(zip(a.seq[preadi:preadi+bases], a.query_qualities[preadi:preadi+bases])):
            if q<baseq or b not in base2index:
                continue
            block[ii*basesize+base2index[b]+i] += 1
        yield prefi, block

def is_antisense(a):
    """Return 1 if read pair is from antisense strand"""
    #if a.is_read1 and a.is_reverse or a.is_read2 and not a.is_reverse or not a.is_paired and a.is_reverse:
    if a.is_reverse:
        if a.is_read1 or not a.is_paired:
            return 1
    elif a.is_read2:
        return 1
    return 0
        
def bam2calls(bam, ref, start, end, mapq=15, baseq=20, offset=33):
    """Return 2D array of basecalls from BAM file, as follows
    - 1D positions from start to end
    - 2D base counts for ACGT from sense and antisense strand at given position
    """
    sam = pysam.AlignmentFile(bam)
    # ACGT x2 for each strand
    basesize = 2*len(alphabet)
    n =  basesize * (end-start+1)
    calls = np.zeros(n, dtype="int64") 
    # stop if ref not in sam file
    if ref not in sam.references:
        return calls.reshape((end-start+1, basesize))
    pa = None  
    for a in sam.fetch(ref, start, end):
        if is_qcfail(a, mapq) or is_duplicate(a, pa): continue
        pa = a
        # get transcript strand
        i = 0 # for +/for i == 0; for -/rev i==len(alphabet)+1
        if is_antisense(a):
            i = len(alphabet)
        for refi, block in get_blocks(a, start, end, baseq, i, basesize):
            s, e = basesize*(refi-start), basesize*(refi-start)+len(block)
            calls[s:e] += block
    return calls.reshape((end-start+1, basesize))

def get_combined_calls(bams, ref, start, end, mapq, baseq, stranded=0):
    """Combine basecalls from several files"""
    parsers = (bam2calls(bam, ref, start, end, mapq, baseq) for bam in bams)
    for call in np.sum(parsers, axis=0): 
        if stranded:
            yield (call[:len(alphabet)], call[len(alphabet):])
        else:
            yield (call[:len(alphabet)] + call[len(alphabet):],)

def fasta2calls(fastafn, ref, start, end, cov=100):
    """Return list of basecalls from FastA file."""
    fasta = pysam.FastaFile(fastafn)
    if ref not in fasta.references:
        raise StopIteration
    for b in fasta.fetch(ref, start, end):
        call = [0]*len(alphabet)
        if b in base2index:
            call[base2index[b]] += cov
        yield call

def get_calls(dna, fasta, stranded, ref, start, end, mapq, baseq, minDepth):
    """Return basecalls from multiple BAM (& FastA) file(s)"""
    # define strands
    if not stranded:
        strands = "."  # unstranded
    elif stranded=="firststrand":
        strands = "+-" # dUTP, NSR, NNSR
    else:
        strands = "-+" # Illumina or Standard Solid
    # get parsers    
    refparser = fasta2calls(fasta, ref, start, end)
    dnaparser = get_combined_calls(dna, ref, start, end, mapq, baseq, stranded)
    # process
    for pos, (refcall, dnacalls) in enumerate(zip(refparser, dnaparser), start+1):
        if sum(refcall)<minDepth:
            continue
        #print ref, pos, refcall, dnacalls
        strand_info = sorted(zip(dnacalls, strands), key=lambda x: sum(x[0]), reverse=1)
        yield ref, pos, refcall, strand_info

def get_allele_freqs(counts, minFreq=0.03, minCount=3):
    """Return two lists: alleles that passed filtering and their frequencies"""
    bases, freqs = [], []
    for c, b in zip(counts, alphabet):
        freq = 1.*c/sum(counts)
        # skip if alt base calling if less than 3 reads or low freq
        if c >= minCount and freq >= minFreq:
            bases.append(b)
            freqs.append(freq)
    return bases, freqs        
        
def bam2heterozygous(dna, fasta, stranded, minDepth, minDNAfreq, mapq, baseq, verbose, ref, start, end):
    """Return RNA editing positions"""
    if verbose:
        region = "%s:%s-%s"%(ref, start, end)
        sys.stderr.write(" %s    \r"%region)
    parser = get_calls(dna, fasta, stranded, ref, start, end, mapq, baseq, minDepth)
    for contig, pos, refbases, strand_info in parser:
        baseRef, refFreq = get_allele_freqs(refbases, minDNAfreq)
        if not baseRef: continue
        
        for sbases, strand in strand_info:
            if sum(sbases)<minDepth: continue
            
            bases, freqs = get_allele_freqs(sbases, minDNAfreq)
            if not bases or bases==baseRef: continue

            fb = sorted(((f,b) for b, f in zip(bases, freqs)), reverse=1)
            
            alts = "/".join(b for f, b in fb)
            freqs = "/".join("%.3f"%f for f, b in fb)
            #print ref, pos, baseRef, sbases
            yield contig, pos, sum(sbases), alts, freqs

def init_args(*args):
    global dna, fasta, stranded, minDepth, minDNAfreq, mapq, bcq, verbose
    dna, fasta, stranded, minDepth, minDNAfreq, mapq, bcq, verbose = args
    
def worker(args):
    global dna, fasta, stranded, minDepth, minDNAfreq, mapq, bcq, verbose
    ref, start, end = args
    totdata = []
    for data in bam2heterozygous(dna, fasta, stranded, minDepth, minDNAfreq,
                               mapq, bcq, verbose, ref, start, end):
        totdata.append(data)
    return totdata
            
def get_regions(bams, step=1000000):
    """Return chromosome regions covered by at least mincov"""
    sam = pysam.Samfile(bams[0])
    references, lengths = sam.references, sam.lengths
    for ref, length in zip(references, lengths):
        for s in range(0, length, step):
            yield ref, s, s+step
            
def logger(info, add_timestamp=1, add_memory=1, out=sys.stderr):
    """Report nicely formatted stream to stderr"""
    memory = timestamp = ""
    if add_timestamp:
        timestamp = "[%s] "%datetime.ctime(datetime.now())
    if add_memory:
        selfmem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.
        childrenmem = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.
        memory = " [memory: %7.1f Mb]"%(childrenmem + selfmem, ) #; %7.1f Mb self
    out.write("%s%s%s\n"%(timestamp, info, memory))
            
def main():
    import argparse
    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.15b')
    parser.add_argument("-o", "--output", required=True,  help="output file")
    parser.add_argument("-d", "--dna", nargs="*", default = [],  help="input DNA-Seq BAM file(s)")
    parser.add_argument("-f", "--fasta", default = None,  help="reference FASTA file")
    parser.add_argument("-s", "--stranded", "-fr-secondstrand", default=False, action="store_true", 
                        help="stranded RNAseq libraries ie. Illumina or Standard Solid")
    parser.add_argument("-fr-firststrand", default=False, action="store_true", 
                        help="stranded RNAseq libraries ie. dUTP, NSR, NNSR")
    parser.add_argument("--minDepth", default=5,  type=int,
                        help="minimal depth of coverage [%(default)s]")
    parser.add_argument("--minDNAfreq",  default=0.05, type=float,
                        help="min frequency for genomic base [%(default)s]")
    parser.add_argument("-m", "--mapq", default=15, type=int, help="mapping quality [%(default)s]")
    parser.add_argument("--bcq", default=20, type=int, help="basecall quality [%(default)s]")
    parser.add_argument("-t", "--threads", default=1, type=int, help="number of cores to use [%(default)s]")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
     
    # mark stranded protocol
    if o.fr_firststrand:
        o.stranded = "firststrand"
    
    # check if all input files exists
    for fn in o.dna:
        if not os.path.isfile(fn):
            sys.stderr.write("No such file: %s\n"%fn)
            sys.exit(1)

    # check if outfile exists and not empty
    if o.output=="-":
        output = sys.stdout
    elif os.path.exists(o.output) and open(o.output).readline():
        sys.stderr.write("The output file %s exists!\n"%o.output)
        sys.exit(1)
    else:
        output = open(o.output, "w")

    runinfo = " ".join(sys.argv)
    header = "## %s\n# chr\tpos\tcov\tbases\tfreqs\n"%runinfo
    output.write(header)
    output.flush()
    
    logger("Indexing bam file(s)...")
    for fn in o.dna:
        if not os.path.isfile(fn+".bai"):
            cmd = "samtools index %s"%fn
            if o.verbose:
                sys.stderr.write(" %s\n"%cmd)
            os.system(cmd)    

    logger("Genotyping...")
    info = "%s\t%s\t%s\t%s\t%s\n"
    regions = get_regions(o.dna)
    if o.threads<2: # this is useful for debugging
        for ref, start, end in regions:
            #print ref, start, end
            #if start>10000: break
            parser = bam2heterozygous(o.dna, o.fasta, o.stranded, o.minDepth, o.minDNAfreq, \
                                      o.mapq, o.bcq, o.verbose, ref, start, end)
            for data in parser:
                output.write(info%data)
    else:
        initargs = (o.dna, o.fasta, o.stranded, o.minDepth, o.minDNAfreq, o.mapq, o.bcq, o.verbose)
        p = Pool(o.threads, initializer=init_args, initargs=initargs)
        parser = p.imap_unordered(worker, regions)
        for data in parser:
            output.write("".join(info%d for d in data))

    logger("Done!")
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
            