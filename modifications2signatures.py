#!/usr/bin/env python
desc="""Generate RT signatures for modifications
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Warsaw, 1/12/2017
"""

import os, sys
reload(sys)  
sys.setdefaultencoding('utf8')

import re, subprocess
from datetime import datetime
from collections import Counter

from modifications2rna import fasta_parser, table2modifications
import numpy as np

import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
import matplotlib.pyplot as plt

import urllib, urllib2

#find stats of the reads in mpileup
##http://samtools.sourceforge.net/pileup.shtml
read_start_pat = re.compile('\^.')
indel_pat = re.compile('[+-]\d+')

def load_modifications(rna, wt=set('ACGU'), log=sys.stderr):
    """Return dictionary with modifications for each ref.
    Coordinates are 1-based / GTF style"""
    log.write("Loading modifications...\n")
    # parse fasta
    ii = 0
    mods, unmods = {}, {}
    for name, seq in fasta_parser(rna):
        # store bases
        for i, b in enumerate(seq, 1):
            ii += 1
            if b=="_":
                pass
            elif b in wt:
                if b not in unmods:
                    unmods[b] = []
            else:
                if b not in mods:
                    mods[b] = []
                mods[b].append("%s:%s"%(name,i))
    log.write(" %s bases with %s modifications (%s unique)\n"%(ii, sum(map(len, mods.itervalues())), len(mods)))
    return mods, unmods
    
def _remove_indels(alts):
    """Return mpileup without indels and read start/end marks and number of insertions and deletions at given position
    .$....,,,,....,.,,..,,.,.,,,,,,,....,.,...,.,.,....,,,........,.A.,...,,......^0.^+.^$.^0.^8.^F.^].^],
    ........,.-25ATCTGGTGGTTGGGATGTTGCCGCT..
    """
    ends = alts.count('$')
    # But first strip start/end read info.
    starts = read_start_pat.split(alts)
    alts = "".join(starts).replace('$', '')
    ends += len(starts)-1
    # count indels
    indels = {"+": 0, "-": alts.count('*')}
    # remove indels info
    m = indel_pat.search(alts)
    while m:
        # remove indel
        pos = m.end() + int(m.group()[1:])
        # count insertions and deletions
        indels[m.group()[0]] += 1
        alts = alts[:m.start()] + alts[pos:]
        # get next match
        m = indel_pat.search(alts, m.start())
    return alts, indels["+"], indels["-"], ends

def genotype_region(bams, dna, region, minDepth, mpileup_opts, alphabet='ACGT'):
    """Start mpileup"""    
    # open subprocess #'-f', dna, 
    args = ['samtools', 'mpileup', '-r', region] + mpileup_opts.split() + bams
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) #, bufsize=65536)
    # process lines
    data, ends = [], []
    for line in proc.stdout:
        data.append([])
        ends.append([])
        lineTuple = line.strip().split('\t')
        # get coordinate
        contig, pos, baseRef = lineTuple[:3]
        baseRef, refFreq = [baseRef], [1.0]
        samplesData = lineTuple[3:]
        # process bam files
        for alg in samplesData[1::3]:
            # remove indels & get base counts
            alts, ins, dels, _ends = _remove_indels(alg)
            counts = [alts.upper().count(b) for b in alphabet] + [ins, dels]
            data[-1].append(counts)
            ends[-1].append(_ends)
            '''if sum(counts)>=minDepth:
                data[-1].append(counts)
            else:
                data[-1].append([0]*len(counts))'''
    return data, ends

def array2plot(outfn, mod, title, cm, pos, window, width=0.75, alphabet='ACGT+-',
               colors=['green', 'blue', 'orange', 'red', 'grey', 'black']):
    """Genotype positions"""
    fig = plt.figure(figsize=(7, 4+3*len(pos)))
    fig.suptitle(title, fontsize=20)
    ind = np.arange(-window-width/2, window)
    for j in range(cm.shape[0]):
        ax = fig.add_subplot(len(pos), 1, j+1)
        ax.set_title(pos[j])
        ax.set_ylim(0,1)
        # plot stacked barplot
        bottom = np.zeros(len(ind))
        for i in range(len(ind)):
            p = ax.bar(ind, cm[j,:,i], width, label=alphabet[i], color=colors[i], bottom=bottom)
            bottom += cm[j,:,i]
            
    #fig.show()#; sys.exit() #format=fformat, 
    fig.savefig(outfn, dpi=100, orientation='landscape', transparent=False)            
        
    if len(pos)>1:
        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title("collapsed signal for %s"%mod)
        ax.set_ylim(0,1)
        # plot combined plot for all positions
        cmm = cm.mean(axis=0)
        # plot stacked barplot
        bottom = np.zeros(len(ind))
        for i in range(len(ind)):
            p = ax.bar(ind, cmm[:,i], width, label=alphabet[i], color=colors[i], bottom=bottom)
            bottom += cmm[:,i]
        fig.savefig(".".join(outfn.split('.')[:-1])+".collapsed."+outfn.split('.')[-1],
                    dpi=100, orientation='landscape', transparent=False)
    # clear
    fig.clear(); del fig
            
def pos2logo(outdir, mod, c, pos, window=2, alphabet='ACGT', ext="svg"):
    """Store logo for each position"""
    # 
    url = 'http://weblogo.threeplusone.com/create.cgi' # "alphabet": auto
    values = {'sequences': '', 'format': ext, 'stack_width': 'medium', 'stack_per_line': '40',
              'alphabet': "alphabet_dna", 'ignore_lower_case': True, 'unit_name': "bits", 'first_index': '1',
              'logo_start': '1', 'logo_end': str(2*window+1), 'composition': "comp_auto", 'percentCG': '',
              'scale_width': True, 'show_errorbars': True, 'logo_title': '', 'logo_label': '',
              'show_xaxis': True, 'xaxis_label': '', 'show_yaxis': True, 'yaxis_label': '',
              'yaxis_scale': 'auto', 'yaxis_tic_interval': '1.0', 'show_ends': True, 'show_fineprint': True,
              'color_scheme': 'color_auto', 'symbols0': '', 'symbols1': '', 'symbols2': '', 'symbols3': '',
              'symbols4': '', 'color0': '', 'color1': '', 'color2': '', 'color3': '', 'color4': ''}
            
    # combine replicas
    csum = c.sum(axis=2)
    for i, p in enumerate(pos):
        freqs = ["P0\tA\tC\tG\tT\n"]
        for j in range(csum.shape[1]):
            freqs.append("%s\t%s\n"%(str(j).zfill(2), "\t".join(map(str, csum[i][j][:len(alphabet)]))))
        # communicate with server and store png
        values["sequences"] = "".join(freqs)
        data = urllib.urlencode(values).encode("utf-8")
        req = urllib2.Request(url, data)
        response = urllib2.urlopen(req); 
        outfn = os.path.join(outdir, "logo.%s.%s.%s"%(mod, p, ext))
        with open(outfn, "wb") as f:
            im = response.read()
            f.write(im)

def modifications2signatures(outdir, bams, dna, rna, table, minDepth, mpileup_opts, verbose, log=sys.stdout, window=2):
    """Generate RT signatures for modifications"""
    mod2base, mod2name = table2modifications(table)
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
        
    # load modifications
    mods, unmods = load_modifications(rna)
    # write header
    log.write("#code\tmodification\toccurencies\tavg cov\tcov std/avg\tA\tC\tG\tT\tins\tdel\n")
    for mod, pos in mods.iteritems(): #list(mods.iteritems())[-10:]:
        data, ends = [], []
        for p in pos:
            ref = p.split(':')[0]
            i = int(p.split(':')[1])
            s, e = i-window, i+window
            if s<1:
                continue
            region = "%s:%s-%s"%(ref, s, e)
            _data, _ends = genotype_region(bams, dna, region, minDepth, mpileup_opts)
            if len(_data)==2*window+1:
                data.append(_data)
                ends.append(_ends)
        if not data:
            continue
        # normalise 0-1 freq
        c, e = np.array(data), np.array(ends)#; print c.shape, e.shape
        if len(c.shape)<4:
            sys.stderr.write("[WARNING] Wrong shape for %s: %s\n"%(mod, c.shape))
            continue
        csum = c.sum(axis=3, dtype='float')
        csum[csum<1] = 1
        cn = 1.*c/csum[:,:,:,np.newaxis]
        # average over all replicas
        cm = cn.mean(axis=2)
        # mean cov & cov var (stdev / mean)
        cov = csum.mean(axis=2).mean(axis=0) 
        covvar = cov.std() / cov.mean()
        # average over all positions
        cmm = cm.mean(axis=0)
        log.write("%s\t%s\t%s\t%.2f\t%.3f\t%s\n"%(mod, mod2name[mod], len(pos), cov[window], covvar, "\t".join("%.3f"%x for x in cmm[window])))
        # plot base freq
        outfn = os.path.join(outdir, "mods.%s.png"%mod2name[mod]) 
        title = "%s [%s] in %s position(s) (%sx)"%(mod2name[mod], mod, len(pos), cov[window])
        array2plot(outfn, mod2name[mod], title, cm, pos, window)
        # store logo
        try:
            pos2logo(outdir, mod2name[mod], c, pos, window)
        except Exception, e:
            sys.stderr.write("[ERROR][pos2logo] %s\n"%str(e))
                
def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.1') 
    parser.add_argument("-o", "--outdir", default="mod2sig", help="output directory [%(default)s]")
    parser.add_argument("-b", "--bam", nargs="+", help="BAM files to process")
    parser.add_argument("-d", "--dna", required=1, help="DNA FastA")
    parser.add_argument("-r", "--rna", required=1, help="RNA FastA")
    parser.add_argument("-t", "--table", default="modifications.txt", help="modification table [%(default)s]" )
    parser.add_argument("-m", "--mpileup_opts", default="-A -q 15 -Q 20", help="options passed to mpileup [%(default)s]")
    parser.add_argument("--minDepth", default=100, type=int, help="minimal depth [%(default)s]")
#    parser.add_argument("-f", "--minFreq", default=0.8, type=float, help="min frequency of alternative base [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write( "Options: %s\n" % str(o) )

    modifications2signatures(o.outdir, o.bam, o.dna, o.rna, o.table, o.minDepth, o.mpileup_opts, o.verbose)

if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
