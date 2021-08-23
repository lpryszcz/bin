#!/usr/bin/env python2
desc="""Filter pairs with at least one read aligned. 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Mizerow, 25/06/2014
"""
import os, sys
import pysam
from datetime import datetime
from Bio.Seq import Seq

def pair2interleaved_fasta(pair):
    """Report interleaved fasta."""
    return "".join(">%s/%s\n%s\n"%(a.qname, i, a.seq) for i, a in enumerate(pair, 1))

def sam2one_aligned(samfn, out, pair2out, verbose):
    """Parse SAM algs and report pairs with at least on read aligned."""
    k = 0
    if samfn=="-":
        sam = pysam.Samfile(samfn, "r")        
    else:
        sam = pysam.Samfile(samfn)
    pair = []
    for i, a in enumerate(sam, 1):
        if verbose and not i%1e5:
            sys.stderr.write(" %s %s\r"%(i, k))
        #skip if both unmapped or if secondary or supplementary alg
        if a.is_unmapped and a.mate_is_unmapped or a.is_secondary or a.flag & 2048:
            continue
        #reverse complement
        if a.is_reverse:
            a.seq, a.qual = str(Seq(a.seq).reverse_complement()), a.qual[::-1]
        #count correct
        k += 1
        pair.append(a)
        #report
        if len(pair)==2:
            out.write(pair2out(pair))
            pair = []
    #report stats
    info = "%s algs processed. %s [%.2f%s] reads reported.\n"
    sys.stderr.write(info%(i, k, k*100.0/i, '%'))

def main():
    import argparse
    usage   = "bwa mem REF read1 read2 | %(prog)s -v" 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b')
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--sam", default="-",      
                        help="input SAM/BAM stream [stdin]")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream [stdout]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
    #define function converting to output format
    pair2out = pair2interleaved_fasta
    #process
    sam2one_aligned(o.sam, o.output, pair2out, o.verbose)

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
