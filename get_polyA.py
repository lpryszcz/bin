#!/usr/bin/env python3
desc="""Estimate polyA tail lengths

Dependencies: h5py mappy ont_fast5_api scipy numpy

"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Cologne, 15/04/2020
"""

import itertools, h5py, mappy, math, os, scipy, sys, numpy as np
from multiprocessing import Pool
from ont_fast5_api.fast5_interface import get_fast5_file
from datetime import datetime
from pathlib import Path
import pandas as pd
from scipy.signal import savgol_filter, find_peaks
from mod_encode import logger, VERSION
#logger = lambda mssg: sys.stderr.write("%s\n"%mssg.rstrip())

def plot_barcode_polyA(sig, s, e, e2, sig_step, sig_step_peaks_pos, sig_step_peaks_neg):
    """ """
    fig, ax = plt.subplots(1, figsize=(20, 5))
    ax1 = ax.twinx()
    ax.plot(sig)
    ax.axvspan(s, e, color="grey", alpha=.3, label="barcode")
    ax1.plot(sig_step)
    ax1.plot(sig_step_peaks_pos, sig_step[sig_step_peaks_pos], "+", color="red")
    ax1.plot(sig_step_peaks_neg, sig_step[sig_step_peaks_neg], "x", color="blue")        
    ax.axvspan(e, e2, color="red", alpha=.3, label="polyA")
    ax.legend()
    #ax.set_xlim(8000, 11000)    

def get_refined_polyA_len(sig, sig_step, peak_width=150, p=0.001, max_iter=3):
    e = sig.shape[0]
    mid = int(round(e/3)); mid2 = 2*mid
    s, p = scipy.stats.ttest_ind(sig[:mid],sig[mid2:])
    while p<0.001 and max_iter:
        max_iter -= 1
        peak_width = int(0.33 * peak_width)
        #print(max_iter, len(sig), peak_width, sig.std(), sig[:mid].std(), sig[mid2:].std(), p)
        sig_step_peaks, _ = find_peaks(sig_step, width=peak_width)
        if len(sig_step_peaks):
            e = sig_step_peaks[0]
            sig, sig_step = sig[:e], sig_step[:e]
            mid = int(round(e/3)); mid2 = 2*mid
            s, p = scipy.stats.ttest_ind(sig[:mid],sig[mid2:])
    #print(max_iter, len(sig), peak_width, sig.std(), sig[:mid].std(), sig[mid2:].std(), p, np.median(sig[:mid]), np.median(sig[mid2:]))
    return len(sig)
 
def get_barcode_start_end(sig, pos_peak_width=500, neg_peak_width=150, max_start=30000, 
                          max_barcode_neg_peaks=3, plot=False):
    "Return barcode start and end computer through convolution"
    sig = sig[:max_start]
    step = np.hstack((np.ones(len(sig)), -1*np.ones(len(sig))))
    # ~35x faster than np.convolve
    sig_step = scipy.signal.fftconvolve(sig-sig.mean(), step, mode='valid') 
    # get barcode boundaries using min and max
    argmax = np.argmax(sig_step)
    # get barcode boundaries using peaks - more sensitive
    sig_step_peaks_pos, _ = find_peaks(sig_step, width=pos_peak_width) # barcode-polyA peak has to be large
    sig_step_peaks_neg, _ = find_peaks(-1*sig_step, width=neg_peak_width) # polyA-end and barcode start can be smaller
    if argmax in sig_step_peaks_pos[1:3]:
        e = argmax
        prev_peak = sig_step_peaks_pos[np.argwhere(sig_step_peaks_pos[:3]==e)[0]-1]
        if not len(sig_step_peaks_neg[sig_step_peaks_neg>prev_peak]): return e, e+1, e+2, sig_step, False
        s = sig_step_peaks_neg[sig_step_peaks_neg>prev_peak][0]
    else:
        s = sig_step_peaks_neg[0]
        if not len(sig_step_peaks_pos[sig_step_peaks_pos>s]): return s, s+1, s+2, sig_step, False
        e = sig_step_peaks_pos[sig_step_peaks_pos>s][0]
    if not len(sig_step_peaks_neg[sig_step_peaks_neg>e]): return s, e, e+1, sig_step, False
    e2 = sig_step_peaks_neg[sig_step_peaks_neg>e][0] 
    # get refined polyA boundaries
    e2 = e + get_refined_polyA_len(sig[e:e2], -1*sig_step[e:e2], neg_peak_width)
    # check if too many neg peaks in barcode
    ok = np.sum(sig_step_peaks_neg<e) <= max_barcode_neg_peaks    
    if plot: plot_barcode_polyA(sig, s, e, e2, sig_step, sig_step_peaks_pos, sig_step_peaks_neg)
    return s, e, e2, sig_step, ok

def get_read_data(f5file, read_id):
    """Return reads signal"""
    read = f5file.get_read(read_id)
    sig = read.get_raw_data(scale=True) #, start=0, end=50967)
    bcgrp = read.get_latest_analysis("Basecall_1D") #Basecall_1D_000
    stride = read.get_analysis_attributes("%s/Summary/basecall_1d_template"%bcgrp)['block_stride']
    #trace = read.get_analysis_dataset(bcgrp, "BaseCalled_template/Trace"); trace
    fq = read.get_analysis_dataset(bcgrp, "BaseCalled_template/Fastq")
    fq_header, fq_seq, fq_spacer, fq_quals = fq.split('\n')[:4] #; fq_seq    
    moves = read.get_analysis_dataset(bcgrp, "BaseCalled_template/Move")
    move_pos = np.argwhere(moves==1).flatten()
    seg_attr = read.get_analysis_attributes("Segmentation_%s/Summary/segmentation"%bcgrp.split("_")[-1])
    first_sample_template = seg_attr['first_sample_template']
    return sig, fq_seq, moves, move_pos, first_sample_template, stride

def get_effective_bases_from_alignment(aligner, seq):
    """Return effective number of bases and sig start position and end"""
    try:
        a = next(aligner.map(seq, buf=None))
    except:
        return len(seq), 0, "-", 0, 0, 0, "."
    # unload alignment
    q_st, q_en, ref, ref_start, ref_end, mapq, strand = a.q_st, a.q_en ,a.ctg, a.r_st, a.r_en, a.mapq, '+' if a.strand == 1 else '-'
    bases_eff = q_en-q_st
    # for RNA sig and move is reverse
    pos_st = len(seq)-q_en if strand=="+" else q_st
    return bases_eff, pos_st, ref, ref_start, ref_end, mapq, strand
        
def process_fast5(fast5, fasta, mapq, quantile=0.05, pos_peak_width=500, neg_peak_width=150):
    """Plot polyA from segmentation"""
    f5file = get_fast5_file(fast5, mode="r")
    aligner = mappy.Aligner(fasta, preset=str('map-ont'), best_n=1)
    # process reads from multi fast5    
    for read_id in f5file.get_read_ids():#[:200]: 
        sig, seq, moves, move_pos, first_sample_template, stride = get_read_data(f5file, read_id)
        # get barcode, polyA and read start # savgol_filter(sig, 51, 7)
        barcode_start, polyA_start, polyA_end, sig_step, ok = get_barcode_start_end(sig, pos_peak_width, neg_peak_width)
        # get polyA length by avg 
        polyA_mean, polyA_stdev = sig[polyA_start:polyA_end].mean(), sig[polyA_start:polyA_end].std()
        # calculate polyA_len and read_end
        polyA_len = polyA_end - polyA_start
        read_end = first_sample_template+move_pos[-1]*stride
        # correct number of bases, read start & alignment details
        bases, read_start_pos, ref, ref_start, ref_end, mapq, strand = get_effective_bases_from_alignment(aligner, seq)
        # update steps
        move_pos = move_pos[read_start_pos:read_start_pos+bases]
        # calculate polyA from avg sig_per_base
        sig_len = (move_pos[-1]-move_pos[0]) * stride
        sig_per_base = sig_len/bases
        polyA_bases = polyA_len/sig_per_base
        # or stripping quantile of largest & smaller steps
        steps = move_pos[1:]-move_pos[:-1]        
        mins, maxs = np.quantile(steps, [quantile, 1-quantile])
        sig_per_base_quant = steps[np.all([steps<=maxs, steps>=mins], axis=0)].mean() * stride
        polyA_bases_quant = polyA_len/sig_per_base_quant
        # find peaks
        if polyA_len>51:
            yhat = savgol_filter(sig[polyA_start:polyA_end], 51, 7)
            peaks, _ = find_peaks(yhat, width=5)
        else:
            peaks = []
        # calculate bases after polyA in the same window
        m_start = move_pos[0] 
        m_end = int(round(m_start+polyA_len/stride))
        bases_after_polyA = moves[m_start:m_end].sum()
        yield (read_id, ok, bases, polyA_bases, polyA_bases_quant, len(peaks), bases_after_polyA, 
               ref, ref_start, ref_end, mapq, strand, 
               polyA_start, polyA_len, sig_len, polyA_mean, polyA_stdev, fast5)
               
def worker(args):
    import warnings
    warnings.filterwarnings('ignore')
    return [out for out in process_fast5(*args)]

def mod_polyA(indirs, fasta, outfn, mapq, threads, recursive,
              min_polyA_len=300, min_first_sample_template=2000):
    """Process Fast5 files from input directories and
    report estimated lenghts of polyA tails
    """
    # prepare df
    columns = ['read_id', 'filter_pass', 'bases', 'polyA_len_from_mean', 'polyA_len_from_corrected_mean', 'polyA_len_from_peaks', 'polyA_len_from, bases_after_polyA', 
               'ref', 'start', 'end', 'mapq', 'strand', #'starts', 'ends', 'gene_id', 'transcript_id', 
               'polyA_sig_start', 'polyA_sig_len', 'read_sig_len', 'polyA_mean', 'polyA_stdev', 'fast5']
    df = pd.DataFrame(columns=columns)
    # define imap, either pool of processes or map
    if threads>1:
        p = Pool(threads)#, maxtasksperchild=1)
        imap = p.imap_unordered
    else:
        imap = map
    # process input directories with Fast5 files
    logger("Processing %s directories...\n"%len(indirs))
    for indir in indirs:
        if recursive:
            fnames = sorted(map(str, Path(indir).rglob('*.fast5')))
        else:
            fnames = sorted(map(str, Path(indir).glob('*.fast5')))
        logger(" %s with %s Fast5 file(s)...\n"%(indir, len(fnames)))
        # multiprocessing here!
        args = [(fast5, fasta, mapq) for fast5 in fnames]
        parser = imap(worker, args)
        for fi, out in enumerate(parser, 1):
            # append & save
            df = df.append(pd.DataFrame(out, columns=columns), ignore_index=True)
            sys.stderr.write(" %s / %s  %s reads  \r"%(fi, len(fnames), df.shape[0]))
            df.to_csv(outfn, sep="\t", index=False)
        logger("  %s reads processed"%df.shape[0])
    # close pool
    if threads>1: p.terminate()
    logger("Done!")
        
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version=VERSION)   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-o", "--output", default="polyA.tsv.gz", help="output filename [%(default)s]")
    parser.add_argument("-i", "--indirs", nargs="+", help="input directory with Fast5 files")
    parser.add_argument("-r", "--recursive", action='store_true', help="recursive processing of input directories [%(default)s]")
    parser.add_argument("-f", "--fasta", required=1, help="reference FASTA file")
    parser.add_argument("-t", "--threads", default=3, type=int, help="number of cores to use [%(default)s]")
    parser.add_argument("-m", "--mapq", default=15, type=int, help="min mapping quality [%(default)s]")

    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))
    
    # create output dir
    if os.path.isfile(o.output):
        logger("File exists: %s"%o.output)
        sys.exit(1)
    outdir = os.path.dirname(o.output)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
        
    # align
    mod_polyA(o.indirs, o.fasta, o.output, o.mapq, o.threads, o.recursive)

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
