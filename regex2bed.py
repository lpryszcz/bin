#!/usr/bin/env python2
desc="""Report BED intervals matching regex expressions.

NOTE: Support for --bothStrands is experimental. It should handle [] well,
but probably won't handle {} in your regex. 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Warsaw, 9/09/2015
"""

import os, re, sys
from datetime import datetime
from Bio import SeqIO

base2rc = {"A": "T", "T": "A", "C": "G", "G": "C",
           "a": "t", "t": "a", "c": "g", "g": "c",
           "[": "]", "]": "[", "{": "}", "}": "{"}

def regex2reverse_complement(regex):
    """Return reverse regex patterns"""
    nregex = []
    for r in regex:
        nr = "".join(reversed([base2rc[c] if c in base2rc else c for c in r]))
        nregex.append(nr)
    return nregex

def rec2matches(out, r, pat, handle, strand="+", score=0, i=0):
    """Match sequences"""
    for i, m in enumerate(pat.finditer(str(r.seq)), 1):
        s, e = m.span()
        out.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(r.id, s, e, m.group(), score, strand))
    return i

def regex2bed(handle, out, regex, bothStrands, ignoreCase, verbose):
    """Match regex in FASTA and report bed."""
    flags = 0
    if ignoreCase:
        flags = re.IGNORECASE
    # get positive and reverse complement patterns
    pRegex = '|'.join(regex)
    pPat = re.compile(pRegex, flags=flags)
    if bothStrands:
        mRegex = '|'.join(regex2reverse_complement(regex))
        mPat = re.compile(mRegex, flags=flags)
    # report matches
    if verbose:
        sys.stderr.write("Scanning chromosomes for %s %s ...\n#chromosome\tlength\t+ matches\t- matches\n"%(pRegex, mRegex))
    i = pCount = mCount = 0
    pCounts, mCounts = [], []
    for i, r in enumerate(SeqIO.parse(handle, 'fasta'), 1):
        #match plus / minus patterns separately
        pCount = rec2matches(out, r, pPat, handle, "+")
        if bothStrands:
            mCount = rec2matches(out, r, mPat, handle, "-")
        pCounts.append(pCount)
        mCounts.append(mCount)
        if verbose:
            sys.stderr.write("%s\t%s\t%s\t%s\n"%(r.id, len(r), pCount, mCount))

    sys.stderr.write("# %s + and %s - matches in %s chromosomes\n"%(sum(pCounts), sum(mCounts), i))

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0a')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--input", type=file, 
                        help="Genome FASTA file")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-r", "--regex", nargs="+", 
                        help="REGEX to match  [%(default)s]")
    parser.add_argument("-s", "--bothStrands", default=False, action='store_true', 
                        help="scan also in reverse complement")
    parser.add_argument("-c", "--ignoreCase", default=False, action="store_true",
                        help="ignore sequence case")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    regex2bed(o.input, o.output, o.regex, o.bothStrands, o.ignoreCase, o.verbose)

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
