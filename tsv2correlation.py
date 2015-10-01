#!/usr/bin/env python
desc="""Calculate correlation over particular columns
"""
epilog="""Author: l.p.pryszcz@gmail.com
Warsaw, 30/09/2015
"""

import os, sys
import numpy as np
from scipy.stats import stats
from datetime import datetime

def tsv2correlations(handle, out, columns, verbose=0):
    """Report summed gene expression from transcripts"""
    if verbose:
        sys.stderr.write("Parsing input...\n")
    data = [[] for c in columns]
    for l in handle:
        # unload data
        lData = l[:-1].split('\t')
        for i, c in enumerate(columns):
            v = lData[c]
            try:
                v = float(v)
            except:
                v = "-"
            data[i].append(v)

    # calculate correlation
    out.write("col1\tcol2\tvalid elements\trho\tpval\n")
    for i in range(len(columns)-1):
        for j in range(i+1, len(columns)):
            scores = [(a,b) for a,b in zip(data[i], data[j]) if a!="-" and b!="-"]
            #scores2 = columns[j]
            rho, pval = stats.spearmanr(scores)
            out.write("%s\t%s\t%s\t%s\t%s\n"%(i, j, len(scores), rho, pval))
            
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0a')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--input", default=sys.stdin, type=file,  
                        help="input stream    [stdin]")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-c", "--columns",   nargs="+", type=int, 
                        help="select columns (0-based)")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    tsv2correlations(o.input, o.output, o.columns, o.verbose)

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
