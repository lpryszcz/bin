#!/usr/bin/env python
desc="""Report chromosome/contig ploidy

TBA:
- some function fitting ploidy with int
- same for alt alleles
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Warsaw, 13/12/2017
"""

import os, sys, pysam, resource
from datetime import datetime
from multiprocessing import Pool
import numpy as np
from scipy import stats
from itertools import izip

alphabet = "ACGT" # i=insertion d=deletion
base2index = {b: i for i, b in enumerate(alphabet)}
for i, b in enumerate(alphabet.lower()):
    base2index[b] = i

# CIGAR operations
"""Op BAM Description +1Q +1R
M 0 alignment match (can be a sequence match or mismatch) yes yes
I 1 insertion to the reference yes no
D 2 deletion from the reference no yes
N 3 skipped region from the reference no yes
S 4 soft clipping (clipped sequences present in SEQ) yes no
H 5 hard clipping (clipped sequences NOT present in SEQ) no no
P 6 padding (silent deletion from padded reference) no no
= 7 sequence match yes yes
X 8 sequence mismatch yes yes
    """
def _match(refi, readi, bases): return refi+bases, readi+bases, True
def _insertion(refi, readi, bases): return refi, readi+bases, False
def _deletion(refi, readi, bases): return refi+bases, readi, False
def _skip(refi, readi, bases): return refi, readi, False
code2function = {0: _match, 7: _match, 8: _match, 1: _insertion, 6: _insertion,
                 2: _deletion, 3: _deletion, 4: _insertion, 5: _skip}

def get_blocks(a, start, end, baseq, basesize):
    """Return tuple of aligned position of query and reference. INDEL aware."""
    readi, refi = 0, a.pos
    for ci, (code, bases) in enumerate(a.cigar):
        prefi, preadi = refi, readi
        refi, readi, data = code2function[code](refi, readi, bases)
        # skip if current before start
        if refi<=start:
            continue
        # typical alignment part
        if data:
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
                block[ii*basesize+base2index[b]] += 1
            yield prefi, block

def is_duplicate(a, pa):
    """Return True if read is duplicate"""
    if pa and a.pos==pa.pos and a.flag==pa.flag and a.cigarstring==pa.cigarstring and a.isize==pa.isize and a.seq==pa.seq:
        return True

def is_qcfail(a, mapq=15):
    """Return True if alignment record fails quality checks"""
    if a.mapq<mapq or a.flag&3840: # or is_heavily_clippped(a): 
        return True

def get_maxfreq():
    """Return freqbins and maxfreq"""
    freqbins = np.arange(.0, 1.01, 0.02)
    maxfreq = np.zeros(len(freqbins), dtype='uint32')            
    return freqbins, maxfreq
        
def bam2cov_freq(bam, region, minDepth, minAltFreq, mapq=15, baseq=20):
    """Return 2 arrays of per position coverage and max base frequency histogram"""
    ref, start, end = region
    sam = pysam.AlignmentFile(bam)
    # ACGT x2 for each strand
    strandsNo = 1
    basesize = strandsNo*len(alphabet)
    n =  basesize * (end-start+1)
    calls = np.zeros(n, dtype="float64") # for compatibility with list
    # maxfreq
    freqbins, maxfreq = get_maxfreq()
    # stop if ref not in sam file
    if ref not in sam.references:
        if ref.startswith('chr') and ref[3:] in sam.references:
            ref = ref[3:]
        elif 'chr%s'%ref in sam.references:
            ref = 'chr%s'%ref
        else:
            return np.array([], dtype='uint16'), maxfreq
    # process alignments
    pa = None  
    for a in sam.fetch(ref, start, end):
        if is_qcfail(a, mapq) or is_duplicate(a, pa):
            continue
        pa = a
        for refi, block in get_blocks(a, start, end, baseq, basesize):
            s, e = basesize*(refi-start), basesize*(refi-start)+len(block)
            calls[s:e] += block
    # reshape
    calls = calls.reshape((end-start+1, len(alphabet)))
    # get coverage
    cov = np.array(calls.sum(axis=1), dtype='uint16')
    # get freq
    freq = 1.*calls.max(axis=1)[cov>=minDepth] / cov[cov>=minDepth]
    for i in np.digitize(freq[freq<1-minAltFreq], freqbins, right=1):
        maxfreq[i] += 1
    return cov, maxfreq

def worker(args):
    # ignore all warnings
    np.seterr(all='ignore')
    bam, region, minDepth, minAltFreq, mapq, bcq = args
    return bam2cov_freq(bam, region, minDepth, minAltFreq, mapq, bcq)

def bam2regions(bam, chrs=[], maxfrac=0.05, step=100000, verbose=0):
    """Return chromosome windows"""
    regions, refs, lens = [], [], []
    sam = pysam.Samfile(bam)
    references, lengths = sam.references, sam.lengths
    for ref, length in izip(references, lengths):
        # skip if chr not selected or too small
        if chrs and ref not in chrs or length<maxfrac*max(lengths):
            if verbose:
                sys.stderr.write(" %s skipped\n"%ref)
            continue
        refs.append(ref)
        lens.append(length)    
        for s in xrange(1, length+1, step):
            regions.append((ref, s, s+step-1))
    return regions, refs, lens
    
def bam2ploidy(bam, minDepth=10, minAltFreq=0.05, mapq=3, bcq=20, threads=4, chrs=[], minfrac=0.05, verbose=1):
    """Get alternative base coverage and frequency for every bam file"""
    # exit if processed
    outfn = "%s.ploidy.tsv"%bam
    if os.path.isfile(outfn) and open(outfn).readline():
        sys.stderr.write(" Outfile exists or not empty: %s\n"%outfn)
        return
    # make sure indexed
    logger(" %s"%bam)
    if not os.path.isfile(bam+".bai"):
        logger(" Indexing %s..."%fn)
        cmd = "samtools index %s"%fn
        if o.verbose:
            sys.stderr.write(" %s\n"%cmd)
        os.system(cmd)
    # get regions
    regions, refs, lens = bam2regions(bam, chrs, minfrac)
    # this is useful for debugging 
    i = 0
    if threads<2:
        import itertools
        p = itertools
    else:
        p = Pool(threads)
    parser = p.imap(worker, ((bam, r, minDepth, minAltFreq, mapq, bcq) for r in regions))
    # process regions
    ref2stats = {}
    pref, covs = "", []
    freqbins, maxfreq = get_maxfreq()
    for i, ((ref, start, end), (_cov, _maxfreq)) in enumerate(izip(regions, parser), 1):
        sys.stderr.write(" %s / %s %s:%s-%s         \r"%(i, len(regions), ref, start, end))
        if ref!=pref:
            if covs: # store covmedia, mean, std and most common freq
                ref2stats[pref] = get_stats(covs, maxfreq); print pref, ref2stats[pref], maxfreq
            # reset
            pref, covs = ref, []
            freqbins, maxfreq = get_maxfreq()
        maxfreq += _maxfreq
        covs.append(_cov)
    # process last output
    if covs:
        ref2stats[pref] = get_stats(covs, maxfreq)
    # get min cov ==> ploidy 1
    covstats = np.array([ref2stats[ref][1] for ref in refs])
    mincov = min(covstats[covstats>covstats.mean()*0.1])
    ploidy = covstats / mincov
    # report
    with open(outfn, "w") as out:
        out.write("#ref\tlen\tcov\tploidy\tmost_common_freq\n")
        for r, l, c, p in izip(refs, lens, covstats, ploidy):
            out.write("%s\t%s\t%.2f\t%.2f\t%s\n"%(r, l, c, p, ",".join(map(str, ref2stats[r][-1]))))

def get_stats(covs, freqs, q=5):
    """Return coverage median, mean and stdev"""
    cov = np.concatenate(covs, axis=0)
    if cov.sum()<100: return 0, 0, 0, []
    # get rid of left / right 5 percentile
    mincov, maxcov = stats.scoreatpercentile(cov, q), stats.scoreatpercentile(cov, 100-q)
    cov = cov[np.all(np.array([cov<maxcov, cov>mincov]), axis=0)]
    if cov.sum()<100: return 0, 0, 0, []
    # most common freq
    mostcommon = [x*2 for x in reversed(freqs.argsort()) if freqs[x]>=0.66*freqs.max()]
    return np.median(cov), cov.mean(), cov.std(), mostcommon
        
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
    parser.add_argument('--version', action='version', version='1.20a')
    parser.add_argument("-i", "--bams", nargs="+", help="input BAM file(s)")
    parser.add_argument("-q", "--mapq", default=3, type=int, help="mapping quality [%(default)s]")
    parser.add_argument("-Q", "--bcq", default=20, type=int, help="basecall quality [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int, help="number of cores to use [%(default)s]")
    parser.add_argument("-c", "--chrs", nargs="*", default=[], help="analyse selected chromosomes [all]")
    parser.add_argument("--minDepth", default=10, type=int,  help="minimal depth of coverage for genotyping [%(default)s]")
    parser.add_argument("--minAltFreq", default=0.1, type=float, help="min frequency for DNA base [%(default)s]")
    parser.add_argument("--minfrac", default=0.05, type=float, help="min length of chr/contig as fraction of the longest chr [%(default)s]")

    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    # check if all files exists
    for fn in o.bams:
        if not os.path.isfile(fn):
            sys.stderr.write("No such file: %s\n"%fn)
            sys.exit(1)

    logger("Processing %s BAM file(s)..."%len(o.bams))
    for bam in o.bams:
        bam2ploidy(bam, o.minDepth, o.minAltFreq, o.mapq, o.bcq, o.threads, o.chrs, o.minfrac, o.verbose)

    logger("Finished!")
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
