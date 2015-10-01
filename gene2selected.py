#!/usr/bin/env python
desc="""Select rows from given files based on first column entry. Ordered as example set. 
"""
epilog="""Author: l.p.pryszcz@gmail.com
Warsaw, 30/09/2015
"""

import os, sys
import numpy as np
from datetime import datetime

def gene2selected(handle, out, examples, header=0, verbose=0):
    """Report summed gene expression from transcripts"""
    if verbose:
        sys.stderr.write("Parsing example IDs...\n")
    rid2id = { l.strip(): i for i, l in enumerate(examples) if l.strip() }
    if verbose:
        sys.stderr.write(" %s loaded.\n"%len(rid2id))
    data = ['\n']*len(rid2id)

    if verbose:
        sys.stderr.write("Parsing input...\n")
    for i, l in enumerate(handle):
        # write header
        if i<header:
            out.write(l)
            continue
        # unload data
        lData = l[:-1].split('\t')
        rid = lData[0]
        if rid in rid2id:
            data[rid2id[rid]] = l
    
    out.write("".join(data[header:]))
            
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
    parser.add_argument("-s", "--selected",   required=True, type=file, 
                        help="annotation gtf")
    parser.add_argument("--header", default=0, type=int, 
                        help="header lines [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    gene2selected(o.input, o.output, o.selected, o.header, o.verbose)

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
