#!/usr/bin/env python
desc="""Fold sequences from BED file and report visual representation. 
 
CHANGELOG:
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw, 1/09/2015
"""

import argparse, os, re, sys, pysam
from datetime import datetime
import subprocess
import numpy as np
# https://github.com/pkerpedjiev/forgi or sudo easy_install -U forgi
#import forgi.graph.bulge_graph as cgb

def init_RNALfold(bin="", verbose=0):
    #open subprocess
    args = ['%sRNALfold'%bin, ]
    if verbose:
        sys.stderr.write("Running %s ...\n" % " ".join(args))
    proc = subprocess.Popen(args, shell=1, stdout=subprocess.PIPE, \
                            stdin=subprocess.PIPE, bufsize=1)
    return proc

def bed2seq(bed, fasta, window):
    """Sequence generator from BED file.
    Return sequence of overlapping BED intervals
    and bed intervals contained withing this sequence.
    """
    pchrom, pe = "", 0
    beds = []
    for l in bed:
        chrom, s, e, name = l[:-1].split('\t')[:4]
        s, e = int(s), int(e)
        if chrom != pchrom or s>pe+window:
            # report
            if beds:
                #print pchrom, pe, chrom, s, e
                # get seq
                _chrom = beds[0][0]
                _s = beds[0][1] - window
                if _s<0:
                    _s=0
                _e = beds[-1][2] + window
                seq = fasta.fetch(_chrom, _s, _e)
                # yield
                yield seq, beds, _chrom, _s, _e
            # reset
            beds = []
        # store previous chrom and end
        beds.append((chrom, s, e, name))
        pchrom, pe = chrom, e
    # report
    if beds:
        # get seq
        _chrom = beds[0][0]
        _s = beds[0][1] - window
        if _s<0:
            _s=0
        _e = beds[-1][2] + window
        seq = fasta.fetch(_chrom, _s, _e)
        # yield
        yield seq, beds, _chrom, _s, _e

def seq2elementary_structure(seq, proc):
    """Return elementary structure representation.
    
    (((((((((...((((((.........))))))........((((((.......))))))..)))))))))
    sssssssssmmmsssssshhhhhhhhhssssssmmmmmmmmsssssshhhhhhhssssssmmsssssssss

    's' indicates a stem pair,
    'm' is a multiloop,
    'h' is a hairpin,
    'i' is an interior loop

    For more details have a look in https://www.biostars.org/p/4300/#71309 
    """
    # get local structures
    #stdout, strerr = proc.communicate(seq)
    #out = stdout.split('\n')[:-3]
    proc.stdin.write(seq+"\n")
    proc.stdin.flush()
    out = [] #proc.stdout.readline().split('\n')[:-3]
    l = ""
    while not l.startswith(" "):
        if l:
            out.append(l)
        l = proc.stdout.readline()
    out = out[:-1]
    # init structure holder
    data = np.zeros((len(out), len(seq)), dtype='S1')
    # add local structures
    for i, l in enumerate(out):
        lData = l.split()
        brackets = lData[0]
        s = int(lData[-1])
        # energy in kcal/mol
        energy = float(lData[-2].strip('()')) 
        # update data
        data[i][s-1:s-1+len(brackets)] = list(brackets)
        '''
        # get elementary representation
        bg = cgb.BulgeGraph()
        bg.from_dotbracket(brackets)
        elementary = bg.to_element_string()
        # update data
        data[i][s-1:s-1+len(elementary)] = list(elementary)
        #'''

    # collapse
    signals = []    
    for c in data.T:
        cFiltered = filter(lambda x: x, c)
        if not cFiltered:
            sig = "-"
        else:
            # get probability of base being ds
            dsProb = np.mean([0 if x in "." else 1 for x in cFiltered])
            sig = "s" if dsProb>=0.5 else "l"
        signals.append(sig)
    return "".join(signals)
        
def bed2rnafold(out, bed, fasta, window, verbose, ViennaPath):
    """ """
    out.write("###\n# %s\n"%bed.name)
    # init RNALfold
    proc = init_RNALfold(bin=ViennaPath)

    # iterate bed entries and get sequences
    Es = Ls = Ss = 0
    for seq, beds, chrom, start, end in bed2seq(bed, fasta, window):
        #
        structure = seq2elementary_structure(seq, proc)
        # mark edited sites
        ls = ss = 0
        editing   = [" "]*len(seq)
        for chrom, s, e, name in beds:
            i = s-start
            editing[i] = "|"
            if   structure[i] == "l":
                ls += 1
            elif structure[i] == "s":
                ss += 1
        editing   = "".join(editing)
        # report
        region = ">%s:%s-%s"%(chrom, start, end)
        header = "%s\t%s\t%s\t%s" % (region, len(beds), ls, ss)
        out.write("%s\n%s\n%s\n%s\n" % (region, seq, structure, editing))
        # stats
        Es += len(beds)
        Ls += ls
        Ss += ss
    # report stats
    sys.stderr.write("%s\t%s\t%s\t%s\n"%(bed.name, Es, Ls, Ss))   
    
def main():

    usage  = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.1')
    parser.add_argument("-b", "--bed", type=file, nargs="+", 
                        help="input BED file(s)")
    parser.add_argument("-f", "--fasta", required=1,  
                        help="reference genome fasta")
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("w"), 
                        help="output stream [stdout]")
    parser.add_argument("-w", "--window", default=25,  type=int,
                        help="window size [%(default)s]")
    parser.add_argument("--ViennaPath", default="", #~/src/ViennaRNA/Progs
                        help="path to Vienna Package [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    # open fasta file
    fasta = pysam.Fastafile(o.fasta)

    # iter bed files
    sys.stderr.write("#name\tediting\tloop\tstem\n")
    for bed in o.bed:
        bed2rnafold(o.out, bed, fasta, o.window, o.verbose, o.ViennaPath)
    
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
