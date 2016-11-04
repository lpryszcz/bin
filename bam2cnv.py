#!/usr/bin/env python
desc="""Report scaffolds by joining contigs based on contact matrix from BAM file. 

TBD:
- allow for fitting contig into gap within large contig?
- distance estimation?
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Bratislava, 4/11/2016
"""

import os, resource, subprocess, sys
import scipy.cluster.hierarchy as sch
import numpy as np
from datetime import datetime
from collections import Counter
from FastaIndex import FastaIndex

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

def _get_samtools_proc(bam, mapq=0, regions=[], skipFlag=3844):
    """Return samtools subprocess"""
    # skip unmapped, secondary, QC fails and supplementary algs
    args = map(str, ['samtools', 'view', "-q", mapq, "-F", skipFlag, bam])
    # add regions
    if regions:
        args += regions
    # start subprocess
    proc = subprocess.Popen(args, stdout=subprocess.PIPE)
    return proc

def get_window(chrom, start, length, windowSize, chr2window):
    """Return window alignment belongs to. """
    end = start + length
    # get middle position
    pos = int(start+round((end-start)/2))
    # get window
    window = pos / windowSize + chr2window[chrom]
    return window
    
def bam2array(windows, windowSize, chr2window, bam, mapq, regions=[], verbose=1):
    """Return contact matrix based on BAM"""
    i = 1
    arrays = []
    if verbose:
        logger("Parsing %s BAM files..."%len(bam))
    for _bam in bam:
        # init empty array
        _arrays = [np.zeros(len(w), dtype='uint32') for w in windows]
        # init samtools
        if verbose:
            logger(" %s"%_bam)
        proc = _get_samtools_proc(_bam, mapq, regions)
        # process reads from given library
        for i, line in enumerate(proc.stdout, i):
            if verbose and not i%1e5:
                sys.stderr.write(" %s \r"%i)
            # unload data
            ref1, start1, mapq, cigar, ref2, start2, insertSize, seq = line.split('\t')[2:10]
            start1, start2, seqlen = int(start1), int(start2), len(seq)
            # update ref if alignment in the same chrom
            #if ref2 == "=": ref2 = ref1
            # skip if contig not present in array
            if ref1 not in chr2window[0]: # or ref2 not in chr2window[0]:
                continue
            # get windows
            for ii in range(len(windowSize)):
                w = get_window(ref1, start1, seqlen, windowSize[ii], chr2window[ii])
                # update contact array
                _arrays[ii][w] += 1
        # add arrays
        arrays.append(_arrays)
        # stop subprocess
        proc.terminate()
    if verbose:
        logger(" %s alignments parsed"%i)
    return arrays

def bam2sv_ref(bam, fasta, out, windowSize, mapq=10, regions=[], verbose=1, ploidy=2, mindiff=0.75):
    """Report likely duplications & deletions from coverage"""
    # generate windows
    windowSize, windows, chr2window = fasta2windows(fasta, windowSize, verbose)

    if os.path.isfile("sv.npz"):
        npy=np.load('sv.npz')
        arrays = npy[npy.files[0]]
    else:    
        arrays = bam2array(windows, windowSize, chr2window, bam, mapq, regions=regions, verbose=verbose)
        # save
        np.savez_compressed(out, arrays)
    
    # filter windows with gaps
    
    # report
    for ii, _windowSize in enumerate(windowSize):
        # get normalised ref coverage
        ref = ploidy*arrays[0][ii]/arrays[0][ii].mean()
        for jj in enumerate(bam[1:], 1):
            s = 2*arrays[jj][ii]/arrays[jj][ii].mean()
            # putative dup
            dups = filter(lambda x: len(x)>=5, consecutive(np.argwhere(s-ref > mindiff).T[0], stepsize=3))
            # putative del
            dels = filter(lambda x: len(x)>=5, consecutive(np.argwhere(ref-s > mindiff).T[0], stepsize=3))
            # get only consecutive 3 and cluster

def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) > stepsize)[0]+1)
        
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
    parser.add_argument("-q", "--mapq", default=10, type=int, 
                        help="min mapping quality for variants [%(default)s]")
    parser.add_argument("-w", "--windowSize", nargs="+", default=[200], type=int,
                        help="window size [%(default)s]")
    parser.add_argument("-r", "--regions", nargs="+", default=[], 
                        help="work only in these genomic regions [all]")
    parser.add_argument("-x", "--bed", default='', 
                        help="BED file with intervals to exclude")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    # process
    bam2sv_ref(o.bam, o.fasta, o.out, o.windowSize, o.mapq, o.regions, o.verbose)

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
    