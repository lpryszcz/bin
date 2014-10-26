#!/usr/bin/env python
desc="""Fetch proteomes of taxa from given group form metaPhOrs. 
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 18/07/2014
"""

import argparse, os, sys
import numpy as np
from datetime import datetime
from Bio import SeqIO

def blastout2identity(out, query, blastfiles, evalue, query_overlap, verbose):
    """Return profile from BLAST output"""
    #get query ids
    queries = []
    q2i, q2len = {}, {}
    for i, r in enumerate(SeqIO.parse(query, 'fasta')):
        q = r.id
        queries.append(q)
        q2i[q] = i
        q2len[q] = len(r)

    #write header
    out.write("#query\ttarget\tidentity\toverlap\n")
    #parse results
    for j, handle in enumerate(blastfiles):
        q2identinty = {}
        for l in handle:
            #orf19.10        M!03795532_GIBZA|FG01284.1|I1RCH0       76.19   \
            #42      10      0       361     402     369     410     1.8e-10 64.0
            q, t, ident, matches, mmatches, indels, qs, qe, ts, te, e, bscore = l.split('\t')
            qs, qe, ident, e = int(qs), int(qe), float(ident), float(e)
            #filter by e-value or query overlap
            if e > evalue:
                continue
            qt = (q, t)
            if qt not in q2identinty:
                q2identinty[qt] = np.zeros(q2len[q])
            #update overlap only if better identity
            if q2identinty[qt][qs-1:qe].min()<ident:
                q2identinty[qt][qs-1:qe] = ident
        #report matches with enough overlap
        for qt in sorted(q2identinty):
            overlap = 1.0*len(q2identinty[qt].nonzero()[0])/len(q2identinty[qt])
            identity = q2identinty[qt][q2identinty[qt].nonzero()].mean()
            if overlap < query_overlap:
                continue
            #report match
            q, t = qt
            out.write("%s\t%s\t%5.2f\t%.3f\n"%(q, t, identity, overlap))

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")
    parser.add_argument('-i', '--blastout', nargs="+", type=file, 
                        help="BLAST output files     [%(default)s]")
    parser.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'),
                        help="output stream          [stdout]")
    parser.add_argument('-f', '--fasta',   type=file, 
                        help="query fasta file       [%(default)s]")
    parser.add_argument('-e', '--evalue',  default=1e-03, type=float,
                        help="E-value cut-off        [%(default)s]")
    parser.add_argument('-q', '--query_overlap', default=0.1, type=float,
                        help="query overlap cut-off  [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    blastout2identity(o.output, o.fasta, o.blastout, o.evalue, o.query_overlap, o.verbose)

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
