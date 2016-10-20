#!/usr/bin/env python
desc="""Fetch proteomes of taxa from given group form metaPhOrs. 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Mizerow, 17/07/2014
"""

import argparse, os, sys
import numpy as np
from datetime import datetime
from Bio import SeqIO

def blastout2profile(out, query, blastfiles, evalue, query_overlap, verbose):
    """Return profile from BLAST output"""
    #get query ids
    queries = []
    q2i, q2len = {}, {}
    for i, r in enumerate(SeqIO.parse(query, 'fasta')):
        q = r.id
        queries.append(q)
        q2i[q] = i
        q2len[q] = len(r)
    
    #prepare empty profile
    profile = np.zeros((len(queries), len(blastfiles)), dtype=int)
    
    #parse results
    for j, handle in enumerate(blastfiles):
        q2overlap = {}
        for l in handle:
            #orf19.10        M!03795532_GIBZA|FG01284.1|I1RCH0       76.19   \
            #42      10      0       361     402     369     410     1.8e-10 64.0
            q, t, ident, matches, mmatches, indels, qs, qe, ts, te, e, bscore = l.split('\t')
            #filter by e-value or query overlap
            if float(e) > evalue:
                continue
            qt = (q, t)
            if qt not in q2overlap:
                q2overlap[qt] = np.zeros(q2len[q], dtype=int)
            #update overlap
            q2overlap[qt][int(qs)-1:int(qe)] = 1
        #report matches with enough overlap
        for qt in q2overlap:
            if q2overlap[qt].mean() < query_overlap:
                continue
            #add match
            q, t = qt
            profile[q2i[q]][j] += 1
            
    #report
    header = "#genes\t" + "\t".join(os.path.basename(f.name) for f in blastfiles)
    out.write(header+"\n")
    for q, p in zip(queries, profile):
        out.write("%s\t%s\n"%(q, "\t".join(map(str, p))))
    #np.savetxt(out, profile, delimiter="\t", header=header, fmt='%i')

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
    parser.add_argument('-e', '--evalue',  default=1e-05, type=float,
                        help="E-value cut-off        [%(default)s]")
    parser.add_argument('-q', '--query_overlap', default=0.3, type=float,
                        help="query overlap cut-off  [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    blastout2profile(o.output, o.fasta, o.blastout, o.evalue, o.query_overlap, o.verbose)

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
