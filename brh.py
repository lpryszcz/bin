#!/usr/bin/env python
desc="""Report orthologs in means of best reciprocal hits (BRH) between two proteomes. 
"""
epilog="""Author: l.p.pryszcz@gmail.com
Bratislava, 11/04/2016
"""

import os, sys
from datetime import datetime
import numpy as np

def psl2best(handle, overlap):
    """Report best match for each query"""
    q2t = {}
    for i, l in enumerate(handle):
        if i<5:
            continue
        ldata = l[:-1].split('\t')
        matches, mismatches, repm, Ns, Qgapc, Qgaps, Tgapc, Tgaps = map(int, ldata[:8])
        strand, q, qsize, qstart, qend, t, tsize, tstart, tend, blocks, bsizes, qstarts, tstarts = ldata[8:]
        qsize, qstart, qend, tsize, tstart, tend = map(int, (qsize, qstart, qend, tsize, tstart, tend))
        # skip if not enough overlap
        if 1.*(qend-qstart)/qsize < overlap or 1.*(tend-tstart)/tsize < overlap:
            continue
        # score or update score if better
        score = matches * 5 + mismatches * -3 + Tgapc * -4 + Tgaps * -1
        if q not in q2t:
            q2t[q] = (t, score)
        elif score > q2t[q][1]:
            q2t[q] = (t, score)
    return q2t

def brh(fastas, out, overlap, verbose):
    """Report best reciprocal hits"""
    fn1, fn2 = fastas
    # run blat
    cmd = "blat -prot %s %s %s"
    out1 = "%s.psl"%fn1
    if not os.path.isfile(out1):
        os.system(cmd%(fn2, fn1, out1))
    out2 = "%s.psl"%fn2
    if not os.path.isfile(out2):
        os.system(cmd%(fn1, fn2, out2))
    # parse psl
    q2t1 = psl2best(open(out1), overlap)
    q2t2 = psl2best(open(out2), overlap)
    
    for q, (t, score) in q2t1.iteritems():
        if q2t2[t] and q==q2t2[t][0]:
            out.write("%s\t%s\n"%(q, t))
        
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0a')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-f", "--fastas", nargs=2, 
                        help="input FASTA files")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("--overlap", default=0.33, type=float, 
                        help="min. overlap between the two proteins [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    brh(o.fastas, o.output, o.overlap, o.verbose)
 
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror)+str(e)+"\n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
