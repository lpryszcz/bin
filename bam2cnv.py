#!/usr/bin/env python
desc="""Report CNVs (duplications & deletions) between reference and sample BAMs. 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Bratislava, 4/11/2016
"""

import os, pysam, resource, subprocess, sys
import numpy as np
from datetime import datetime
from collections import Counter
from FastaIndex import FastaIndex
from multiprocessing import Pool

def logger(message, log=sys.stdout):
    """Log messages"""
    memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    log.write("[%s] %s    [memory: %6i Mb]\n"%(datetime.ctime(datetime.now()), message, memory))

def fasta2windows(fasta, windowSize, verbose=1):
    """Generate windows over chromosomes"""
    # init fasta index
    faidx = FastaIndex(fasta)
    if verbose:
        logger("Parsing FastA file...")
    # generate windows
    windows, chr2window = [[] for w in windowSize], [{} for w in windowSize]
    skipped = []
    for i, c in enumerate(faidx, 1): 
        if i%1e5 == 1:
            sys.stderr.write(' %s   \r'%i)
        # get windows
        size = faidx.id2stats[c][0]
        # skip short contigs    
        if size < min(windowSize):
            skipped.append(size)
            continue
        for i in range(len(windowSize)):
            # skip contig if shorter than given window
            if size < windowSize[i]:
                continue
            # get starting window
            chr2window[i][c] = len(windows[i])
            for start in range(0, size, windowSize[i]):
                windows[i].append((c, start, start+windowSize[i]))
            # skip last entry
            #windows[i][-1] = (c, start, size)
    if verbose:
        logger(' %s bases in %s contigs divided in %s-%s windows. '%(faidx.genomeSize, len(faidx), len(windows[0]), len(windows[-1])))
        if skipped:
            logger('  %s bases in %s contigs skipped.'%(sum(skipped), len(skipped)))
    return windowSize, windows, chr2window

def bam2array(args): 
    bamfn, sizes, windowSize, chr2window, flag = args
    # init empty array
    _arrays = [np.zeros(s, dtype='uint32') for s in sizes]
    # open bam
    sam = pysam.Samfile(bamfn)
    for i, r in enumerate(sam, 1):
        if not i%1e5:
            sys.stderr.write(" %s %s [%.2f%s]      \r"%(bamfn, i, 100.*i/(sam.mapped+sam.unmapped),'%s'))
        # skip unmapped, secondary, duplicates and suplementary alignments or 
        if r.flag & flag or r.reference_name not in chr2window[0]:
            continue
        # get windows
        for ii in range(len(windowSize)):
            # get window
            w = (r.pos + r.alen/2) / windowSize[ii] + chr2window[ii][r.reference_name]
            # and update coverage
            _arrays[ii][w] += 1
    return i, _arrays
              
def get_arrays(sizes, windowSize, chr2window, bam, mapq, threads, flag=3844, verbose=1):
    """Return coverage in wondows based on alignments from BAM"""
    algs = 0
    arrays = []
    if verbose:
        logger("Parsing %s BAM files using %s threads..."%(len(bam), threads))
    iterator = [(bamfn, sizes, windowSize, chr2window, flag) for bamfn in bam]
    p = Pool(threads)
    for i, _arrays in p.imap(bam2array, iterator):
        # add arrays
        arrays.append(_arrays)
        algs += i
    if verbose:
        logger("  %s alignments parsed"%algs)
    return arrays

def get_consecutive(data, stepsize=1):
    """Return consecutive windows allowing given max. step size"""
    return np.split(data, np.where(np.diff(data) > stepsize)[0]+1)
    
def get_cnvs(name, diff, windows, chrs, minlen, mindiff, stepsize, verbose):
    """Return collapsed list of CNVs"""
    # note, a CNV may span across chromosomes, thus need to process chromosomes separately
    cnvs = []
    for s, e in chrs:
        # need to add s to arg
        cnvs += filter(lambda x: len(x)>=minlen, get_consecutive(np.argwhere(diff[s:e] > mindiff).T[0]+s, stepsize))
    # collapse
    coords = []
    for i, w in enumerate(cnvs, 1):
        s, e = w[0], w[-1]
        (chrom, start, end), (_chrom, _start, _end) = windows[s], windows[e]
        if chrom != _chrom: print "[ERROR] Wrong chromosome start and end: ", windows[w[0]], windows[w[-1]]
        mean = np.mean(diff[s:e+1])
        coords.append((chrom, start, _end, name%i, mean))
    return coords
    
def bam2sv_ref(bam, fasta, outbase, windowSize, mapq=10, threads=4, ploidy=2, 
               mindiff=0.75, minlen=10, stepsize=3, flag=3844, verbose=1):
    """Report likely duplications & deletions from coverage"""
    # generate windows
    windowSize, windows, chr2window = fasta2windows(fasta, windowSize, verbose)

    # get array
    if os.path.isfile("%s.npz"%outbase):
        npy = np.load('%s.npz'%outbase)
        arrays = npy[npy.files[0]]
    else:
        sizes = [len(w) for w in windows]
        arrays = get_arrays(sizes, windowSize, chr2window, bam, mapq, threads, flag, verbose)
        # save
        np.savez_compressed(outbase, arrays)
    
    # report
    logger("Reporting CNVs...")
    for ii, _windowSize in enumerate(windowSize):
        outfn = "%s.%s.bed"%(outbase, _windowSize)
        print "## %s bp window to %s\n# fname\t# DUPs\t# DELs"%(_windowSize, outfn)
        # get chr boundaries
        starts = sorted(chr2window[ii].values())+[len(windows[ii])]
        chrs = [(starts[i], starts[i+1]) for i in range(len(starts)-1)]
        #
        cnvs = []
        # get normalised ref coverage
        ref = ploidy * arrays[0][ii] / arrays[0][ii].mean()
        for jj, bamfn in enumerate(bam[1:], 1):
            s = ploidy * arrays[jj][ii] / arrays[jj][ii].mean()
            # putative dups and dels
            dups = get_cnvs("DUP%s_"+bamfn, s-ref, windows[ii], chrs, minlen, mindiff, stepsize, verbose)
            dels = get_cnvs("DEL%s_"+bamfn, ref-s, windows[ii], chrs, minlen, mindiff, stepsize, verbose) 
            print "%s\t%s\t%s"%(bamfn, len(dups), len(dels))
            cnvs += dels + dups#print dels; print dups
        #report
        with open(outfn, "w") as out:
            out.write("#chrom\tstart\tend\tname\tploidy change\n")
            for chrom, s, e, name, mean in sorted(cnvs):
                out.write("%s\t%s\t%s\t%s\t%.3f\n"%(chrom, s, e, name, mean))

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--bam", nargs="+",
                        help="BAM file(s), first being reference")
    parser.add_argument("-f", "--fasta", type=file,
                        help="Genome FastA file")
    parser.add_argument("-o", "--out", default='cnv',
                        help="output name [%(default)s]")
    parser.add_argument("-p", "--ploidy", default=2, type=int,
                        help="ploidy to report CNVs [%(default)s]")
    parser.add_argument("-q", "--mapq", default=10, type=int, 
                        help="min mapping quality for variants [%(default)s]")
    parser.add_argument("-w", "--windowSize", nargs="+", default=[100], type=int,
                        help="window size [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int,
                        help="no. of processors to use [%(default)s]")
    parser.add_argument("-d", "--mindiff", default=0.75, type=float,
                        help="min difference in coverage to call CNV in given window [%(default)s]")
    parser.add_argument("-l", "--minlen", default=10, type=int,
                        help="min. no. of windows to call CNV [%(default)s]")
    parser.add_argument("-s", "--stepsize", default=3, type=int,
                        help="max distance between consecutive windows within CNV [%(default)s]")
    parser.add_argument("-F", "--flag", default=3844, type=int,
                        help="bitflag of alignments to filter out [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    # process
    bam2sv_ref(o.bam, o.fasta, o.out, o.windowSize, o.mapq, o.threads, o.ploidy, 
               o.mindiff, o.minlen, o.stepsize, o.flag, o.verbose)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
    