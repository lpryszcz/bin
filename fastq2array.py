#!/usr/bin/env python
desc="""Produce contact matrix from SAM/BAM. 

Make sure there are no spaces in your FastA chromosome/contig names, as this will collapse the program.
You can remove this with: 
cut -f1 -d" " genome.fa > genome.corr.fa
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Bratislava, 19/10/2016
"""

import gzip, os, resource, sys
import commands, os, subprocess, sys
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter, MaxNLocator
from FastaIndex import FastaIndex

def logger(message, log=sys.stdout):
    """Log messages"""
    memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    log.write("[%s] %s    [memory: %6i Mb]\n"%(datetime.ctime(datetime.now()), message, memory))

def fasta2windows(fasta, windowSize, verbose):
    """Generate windows over chromosomes"""
    # make sure no spaces in FastA
    fn = fasta.name
    grepout = commands.getoutput('grep ">" %s | grep -m1 " "'%fn)
    if grepout:
        logger("Removing spaces from FastA...")
        os.system('mv %s %s.bck && cut -f1 -d" " %s.bck > %s'%(fn, fn, fn, fn))
    # 
    if verbose:
        logger("Parsing FastA file...")
    # init fasta index
    faidx = FastaIndex(fasta)
    # filter windows so they are smaller than largest chr and withing reasonalbe range toward genome size
    maxchrlen = max(faidx.id2stats[c][0] for c in faidx)
    windowSize = filter(lambda x: 1000*x<maxchrlen and 1000*x<0.01*faidx.genomeSize and 1000*x>0.00001*faidx.genomeSize, windowSize)
    if verbose:
        logger(" selected %s windows [kb]: %s"%(len(windowSize), str(windowSize)))
    windowSize = [w*1000 for w in windowSize]
    # generate windows
    windows, chr2window = [[] for w in windowSize], [{} for w in windowSize]
    genomeSize = 0
    base2chr = {}
    skipped = []
    for i, c in enumerate(faidx, 1): #.sort(minLength=minLength)
        if i%1e5 == 1:
            sys.stderr.write(' %s   \r'%i)
        # get windows
        size = faidx.id2stats[c][0]
        # skip short contigs    
        if size < min(windowSize):
            skipped.append(size)
            continue
        for i in range(len(windowSize)):
            # get starting window
            chr2window[i][c] = len(windows[i])
            for start in range(0, size, windowSize[i]):
                windows[i].append((c, start, start+windowSize[i]))
            # update last entry end
            windows[i][-1] = (c, start, size)
        # get chromosome tick in the middle    
        base2chr[genomeSize+size/2] = c
        # update genomeSize
        genomeSize += size
    if verbose:
        logger(' %s bases in %s contigs divided in %s-%s windows. '%(faidx.genomeSize, len(faidx), len(windows[0]), len(windows[-1])))
        if skipped:
            logger('  %s bases in %s contigs <%sbp skipped.'%(sum(skipped), len(skipped), min(windowSize)))
    return windowSize, windows, chr2window, base2chr, faidx.genomeSize
        
def _get_snap_proc(fn1, fn2, ref, cores, verbose, log=sys.stderr, largemem=0):
    """Return snap-aligner subprocess.
    bufsize: 0 no buffer; 1 buffer one line; -1 set system default.
    """
    # create genome index
    idxfn = ref + ".snap"
    idxcmd = "snap-aligner index %s %s" % (ref, idxfn)
    if largemem:
        idxfn += ".large"
        idxcmd += " -large"
    if not os.path.isdir(idxfn):
        if verbose:
            log.write(" Creating index...\n  %s\n" % idxcmd)
        idxmessage = commands.getoutput(idxcmd)
        log.write(idxmessage)    
    # skip mate rescue
    args = ['snap-aligner', 'paired', idxfn, fn1, fn2, '--b', '-t', str(cores), '-o', '-sam', '-']
    if verbose:
        log.write( "  %s\n" % " ".join(args) )
    #select ids
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=log)
    return proc
    
def _unload_sam(sam):
    return sam[0], int(sam[1]), sam[2], int(sam[3]), int(sam[4]), len(sam[9])

def parse_sam(handle):
    """Return tuple representing entries from SAM."""
    q1 = q2 = ""
    for l in handle:
        l = l.strip()
        if not l or l.startswith('@'):
            continue
        sam = l.split('\t')
        #first in pair
        if int(sam[1]) & 64:
            #skip multiple matches
            if sam[0] == q1:
                continue
            q1, flag1, ref1, start1, mapq1, len1 = _unload_sam(sam)
        else:
            #skip multiple matches
            if sam[0] == q2:
                continue
            q2, flag2, ref2, start2, mapq2, len2 = _unload_sam(sam)
        #report
        if q1 == q2:
            yield q1, flag1, ref1, start1, mapq1, len1, q2, flag2, ref2, start2, mapq2, len2

def get_window(chrom, start, length, flag, windowSize, chr2window):
    """Return window alignment belongs to. """
    if flag & 16:
        end = start
        start += length
    else:
        end = start + length
    # get middle position
    pos = int(start+round((end-start)/2))
    # get window
    window = pos / windowSize + chr2window[chrom]
    return window

def sam2array(a, windowSize, chr2window, outfn, fq1, fq2, ref, cores, mapq, upto, verbose):
    """Process SAM entries and add them to contact array."""
    """Sort contigs starting from the longest and remove too short"""
    #elif 
    # run bwa for all libs
    alglog = open(outfn+".log", "w")
    proc = _get_snap_proc(fq1, fq2, ref, cores, verbose, alglog)
    # process reads from given library
    j = 0
    for i, data in enumerate(parse_sam(proc.stdout), 1):
        q1, flag1, ref1, start1, mapq1, len1, q2, flag2, ref2, start2, mapq2, len2 = data
        if upto and i>upto:
            break
        if verbose and not i%1e4:
            sys.stderr.write(" %s \r"%i)
        # skip low quality alignments
        if mapq:
            if mapq1 < mapq or mapq2 < mapq:
                continue
        j += 1
        if ref1 not in chr2window[0] or ref2 not in chr2window[0]:
            continue
        # get windows
        for ii in range(len(windowSize)):
            w1 = get_window(ref1, start1, len1, flag1, windowSize[ii], chr2window[ii])
            w2 = get_window(ref2, start2, len2, flag2, windowSize[ii], chr2window[ii])
            # matrix is symmetrix, so make sure to add only to one part
            if w2 < w1:
                w1, w2 = w2, w1
            # update contact array
            a[ii][w1][w2] += 1
    if verbose:
        info = " %s-%s:"%(fq1, fq2)
        if i: # %s in different windows [%.2f%s]. k, k*100.0/i, '%'
            info += " %s pairs. %s passed filtering [%.2f%s].\n" % (i, j, j*100.0/i, '%')
        else:
            info += " No pairs were aligned!\n"
        sys.stderr.write(info)
    # stop subprocess
    proc.terminate()
    return a

def plot(outfn, a, genomeSize, base2chr, _windowSize, dpi=300):
    """Save contact plot"""
    
    def format_fn(tick_val, tick_pos):
        """Mark axis ticks with chromosome names"""
        if int(tick_val) in base2chr:
            return base2chr[int(tick_val)]
        else:
            sys.stderr.write("[WARNING] %s not in ticks!\n"%tick_val)
            return ''
            
    # invert base2chr
    base2chr = {genomeSize-b: c for b, c in base2chr.iteritems()}
    # start figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("Contact intensity plot [%sk]"%(_windowSize/1000,))
    # label Y axis with chromosome names
    if len(base2chr)<50:
        ax.yaxis.set_major_formatter(FuncFormatter(format_fn))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        plt.yticks(base2chr.keys())
        ax.set_ylabel("Chromosomes")
    else:
        ax.set_ylabel("Genome position")
    # label axes
    ax.set_xlabel("Genome position")        
    plt.imshow(a+1, cmap=cm.hot, norm=LogNorm(), extent=(0, genomeSize, 0, genomeSize))# 
    plt.colorbar()
    # save
    fig.savefig(outfn+".svg", dpi=dpi, papertype="a4")

def fastq2array(fasta, fastq, outfn, windowSize, mapq=10, cores=1,
                upto=float('inf'), dpi=300, dtype='uint16', verbose=1):
    """Convert SAM to SSPACE TAB file."""
    # get windows
    windowSize, windows, chr2window, base2chr, genomeSize = fasta2windows(fasta, windowSize, verbose)
    if verbose:
        logger("Mapping reads...")
    # init contact matrices
    arrays = [np.zeros((len(w), len(w)), dtype=dtype) for w in windows]
        
    # and populate with reads
    for fq1, fq2 in zip(fastq[0::2], fastq[1::2]):
        arrays = sam2array(arrays, windowSize, chr2window, outfn, fq1, fq2, fasta.name,
                           cores, mapq, upto, verbose)
        
    # save
    if verbose:
        logger("Saving & plotting...")
    for a, _windowSize,_windows in zip(arrays, windowSize, windows):
        _outfn = outfn + ".%sk"%(_windowSize/1000,)
        if verbose:
            logger(" %s"%_outfn)        
        # save windows
        with gzip.open(_outfn+".windows.tab.gz", "w") as out:
            out.write("\n".join("\t".join(map(str, w)) for w in _windows)+"\n")    

        #save normalised
        with open(_outfn+".npz", "w") as out:
            np.savez_compressed(out, a)

        if len(_windows)<1e4:
            plot(_outfn, a, genomeSize, base2chr, _windowSize, dpi)
        elif verbose:
            sys.stderr.write("[WARNING] Very large matrix (%s x %s). Skipped plotting!\n"%(len(_windows), len(_windows)))
    if verbose:
        logger("Done!")

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--fastq", nargs="+", help="FastQ/A files")
    parser.add_argument("-f", "--fasta", type=file, help="Genome FastA file")
    parser.add_argument("-o", "--output", required=1, help="output name")
    parser.add_argument("-w", "--windowSize", nargs="+", default=[10000, 1000, 500, 100, 50, 20, 10, 5, 1], type=int,
                        help="window size in kb [%(default)s]")
    parser.add_argument("-m", "--mapq", default=10, type=int,
                        help="mapping quality [%(default)s]")
    parser.add_argument("-u", "--upto", default=0, type=float,
                        help="process up to this number of reads from each library [all]")
    parser.add_argument("-t", "--threads", default=4, type=int,
                        help="number of processes to use [%(default)s]")
    parser.add_argument("-d", "--dpi", default=300, type=int,
                        help="output images dpi [%(default)s]")
    parser.add_argument("--dtype", default='float32', 
                        help="numpy array data type (try uint16 if MemoryError for 100+k windows) [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
    # check if even number of FastQ files
    if len(o.fastq)%2:
        sys.stderr.write("Number of FastQ/A files has to be even!\n")
        sys.exit(1)
        
    # check if output not present
    fileexists = 0
    for w in o.windowSize:
        outfn = "%s.%sk.npz"%(o.output, w)
        if os.path.isfile(outfn):
            sys.stderr.write("Output file exists: %s !\n"%outfn)
            fileexists = 1
    if fileexists:
        sys.exit(1)

    # create outdir
    if not os.path.isdir(os.path.dirname(o.output)):
        os.makedirs(os.path.dirname(o.output))
        
    # process
    fastq2array(o.fasta, o.fastq, o.output, o.windowSize, o.mapq, o.threads,
                o.upto, o.dpi, o.dtype, o.verbose)

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