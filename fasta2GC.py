#!/usr/bin/env python
desc="""Report GC content.
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 12/11/2013
"""

import argparse, sys
from Bio import SeqIO
from Bio.SeqUtils import GC
     
def main():
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-i", "--input",     default=sys.stdin, type=file, 
                        help="fasta stream    [stdin]")
                         
    o = parser.parse_args()

    for r in SeqIO.parse(o.input, 'fasta'):
        print " %s [%s bp]: %.2f%s" %(r.id, len(r), GC(r.seq), '%')
	
if __name__=='__main__': 
    main()
