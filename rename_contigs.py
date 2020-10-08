#!/usr/bin/env python
desc="""Rename FastA headers with unique names. 

By default contigXXXXXX is used, where XXXXXX is the number of contig within FastA file. 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 8/10/2020
"""

import os, sys
from datetime import datetime

def rename_contigs(name="contig", handle=sys.stdin, out=sys.stdout):
    """Rename FastA headers"""
    i = 0
    for l in handle:
        # replace headers
        if l.startswith(">"):
            i += 1
            if not i%100: 
                sys.stderr.write(" %s \r"%i)
            out.write(">{}{:06d}\n".format(name, i))
        else:
            out.write(l)
        
def main():
    import argparse
    usage = "%(prog)s -v"
    parser = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', action='version', version='1.0b')
    #parser.add_argument("-v", "--verbose", action="store_true", help="verbose")
    parser.add_argument("-i", "--input", default=sys.stdin, type=argparse.FileType('r'),
                        help="input stream [stdin]")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream [stdout]")
    parser.add_argument("-n", "--name", default="contig",
                        help="contig basename [%(default)s]") 

    o = parser.parse_args()
    #if o.verbose: sys.stderr.write("Options: %s\n"%str(o))
        
    rename_contigs(o.name, o.input, o.output)
 
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
