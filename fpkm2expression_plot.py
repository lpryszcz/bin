#!/usr/bin/env python
desc="""Generate expression plot similar to cummeRbund from .sf expression values.

CHANGELOG:
v1.3
- calculate FPKM for salmon output

"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Mizerow, 8/05/2015
"""

import os, sys, time
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

#colors = ['Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd']
colors = ('b','c','r','g','y', 'm','b','c','r','g','y','m', 'b','c','r','g','y','m')

def load_gene2transcripts(transcripts):
    """Return gene2transcripts.
    Supports ensembl fasta files.
    """
    # ENSDART00000007748 ensembl:known chromosome:Zv9:21:822304:832471:-1 gene:ENSDARG00000016476
    gene2transcripts = {}
    for line in transcripts:
        if line.startswith(">"):
            lData = line[1:].split()
            tid = lData[0]
            gid = filter(lambda x: x.startswith("gene:"), lData)[0].split(':')[1]
            if gid not in gene2transcripts:
                gene2transcripts[gid] = [tid]
            else:
                gene2transcripts[gid].append(tid)
    return gene2transcripts

def parse_sf(handles, conditions, transcripts, genes=[]):
    """Parse multiple .sf handles and return expression info
    for all trascripts of each gene."""
    # get gene2transcripts
    gene2transcripts = load_gene2transcripts(transcripts)
    # get tid2gid
    tid2gid = {} #tid: gid for
    for gid, tids in gene2transcripts.iteritems():
        for tid in tids:
            tid2gid[tid] = gid
    # prepare
    gene2fpkms = {gid: [[0]*len(tids) for i in range(len(conditions))] \
                  for gid, tids in gene2transcripts.iteritems()}
    # parse handles
    for i, handle in enumerate(handles):
        for line in handle:
            if line.startswith(('#','Name\t')):
                continue
            # unload line
            lData = line[:-1].split('\t')
            if len(lData)==5:
                tid, length, tpm, fpkm, reads = lData
            # salmon v4.1+ https://github.com/COMBINE-lab/salmon/releases
            else:
                tid, length, tpm, reads = lData
            # get gid
            gid = tid2gid[tid]
            if genes and gid not in genes:
                continue
            # get fpkm
            length, tpm = float(length), float(tpm)
            fpkm = tpm*1000.0 / length
            # store fpkm for given condition and transcript
            gene2fpkms[gid][i][gene2transcripts[gid].index(tid)] = float(fpkm)
    # yield data
    for gid, fpkms in gene2fpkms.iteritems():
        if genes and gid not in genes:
            continue            
        yield gid, gene2transcripts[gid], fpkms
        
def fpkm2expression_plot(genes, handle, out, transcripts, FPKMfrac, verbose):
    """Parse expression and generate expression plots."""
    if not os.path.isdir(out):
        os.makedirs(out)
    # get conditions
    ## salmon/RZE024/quant.sf -> RZE024
    conditions = [h.name.split("/")[-2].split('_')[-1] for h in handle]
    # get parser
    parser = parse_sf(handle, conditions, transcripts, genes)
    # the x locations for the groups
    ind = np.arange(len(conditions))+.25
    # the width of the bars: can also be len(x) sequence
    width = 0.75
    for i, (gid, transcripts, fpkms) in enumerate(parser, 1):
        if verbose:
            sys.stderr.write(" %s %s %s\n"%(i, gid, len(transcripts)))
        #start figure
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # prepare data
        fpkms = np.array(fpkms).transpose()
        cumFPKM = sum(sum(fpkms))
        bottom = np.zeros(len(conditions))
        ji = 0
        for j, _fpkms in enumerate(fpkms):
            if sum(_fpkms) < FPKMfrac*cumFPKM:
                sys.stderr.write("  Skipped %s with %.2f FPKM\n"%(transcripts[j], sum(_fpkms)))
                continue
            if ji >= len(colors):
                sys.stderr.write("[WARNING] Too many transcripts to plot!\n")
                break
            ax.bar(ind, _fpkms, width, bottom=bottom, label=transcripts[j], color=colors[ji])
            plt.xticks(ind+width/2., conditions, rotation=90, fontsize=9)
            #plt.xlabel("bound             unbound             total")
            ax.set_ylabel("FPKM")
            # update handles
            bottom += _fpkms
            ji += 1
        ax.legend(loc="best", fontsize=8)
        ax.set_title(gid)
        #plt.show()
        fig.savefig(os.path.join(out,"%s.svg"%gid), papertype="A3")

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.3')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--fpkm", default=sys.stdin, type=file, nargs="+", 
                        help="isoforms.fpkm_tracking or .sf file(s) [stdin]")
    parser.add_argument("-g", "--genes", nargs="+", 
                        help="genes to plot")
    parser.add_argument("-o", "--output", default="expression_plot",  
                        help="output stream   [%(default)s]")
    parser.add_argument("-t", "--transcripts", default="", type=file, 
                        help="transcripts file; needed to get gene2transcripts for .sf input [%(default)s]")
    parser.add_argument("-f", "--frac", default=0.05, type=float, 
                        help="ignore transcripts with expression below [%(default)s] of gene expression")
   
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    fpkm2expression_plot(o.genes, o.fpkm, o.output, o.transcripts, o.frac, o.verbose)
 
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
