#!/usr/bin/env python
desc="""Report SAM entries with mismatches >= m.

"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 26/10/2012
"""

import argparse, os, re, sys
from Bio      import Seq
from datetime import datetime
import numpy as np

cigarPat = re.compile('\d+\w')
    
def cigar2aligned(cigar, flag, aligned):
    """ """
    # unload cigar
    blocks = cigarPat.findall(cigar)
    # reverse cigar for reverse alignments
    if int(flag) & 16:
        blocks.reverse()
    # update aligned parts of the read
    i = 0
    for c in blocks:
        bases, s = int(c[:-1]), c[-1]
        # N D - affect ref
        # M S H I - affect read
        if s in ('S','H','I'):
            i += bases
        elif s=="M":
            aligned[i:i+bases] = 1
            i += bases
    return aligned

def sam2unaligned(handle, out, algTh, verbose):
    """SAM entries generator.
    Return read sequence and total aligned fraction over all reported alignments."""
    # keep track of aligned fraction and previous reads
    aligned = np.zeros(0, dtype='uint8')
    pqname, pflag, psam = "", "", []
    for sam in handle:
        # skip header
        if sam.startswith("@"):
            continue
        # read sam
        qname, flag, tname, pos, mapq, cigar, mrnm, \
        mpos, tlen, seq, qual = sam.split("\t")[:11]
        # if new read, or mate from the same read
        # report previous and start new read caching
        if qname != pqname or int(flag)&64 != int(pflag)&64:
            # report poorly aligned
            if len(aligned) and aligned.mean() < algTh:
                out.write(psam[0])
                #if len(psam)>2: print '\n> ' + '\n> '.join(psam) + '\n', aligned.mean(), aligned
            # reset aligned - memory optimised
            aligned = np.zeros(len(seq), dtype='uint8')
            # store new values & reset psam
            pqname, pflag = qname, flag
            psam = []
        # store sam
        psam.append(sam)
        # update aligned for given read
        aligned = cigar2aligned(cigar, flag, aligned)
    # report poorly aligned
    if len(aligned) and aligned.mean() < algTh:
        out.write(psam[0])
        
def main():

    usage  = "usage: bwa mem REF fastq | %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", "--input", default=sys.stdin, type=file, 
                        help="output stream")
    parser.add_argument("-o", "--out",     default=sys.stdout, type=argparse.FileType("w"), 
                        help="output stream")
    parser.add_argument("-a", "--algTh",   default=0.33, type=float,
                        help="max. fraction of read aligned  [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    sam2unaligned(o.input, o.out, o.algTh, o.verbose)
        
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
    