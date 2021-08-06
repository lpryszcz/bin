#!/usr/bin/env python3
desc="""Report FastQ from basecalled Fast5 file(s).
Originally from https://github.com/lpryszcz/Pszczyna 

Dependencies: ont_fast5_api
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 30/00/2020
"""

import os, sys
from datetime import datetime
from ont_fast5_api.fast5_interface import get_fast5_file

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='0.10a')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument("-i", "--fast5", nargs="+", help="input Fast5 file(s)")
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("w"), help="output stream [stdout]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    for fn in o.fast5:
        seqs = []
        f5file = get_fast5_file(fn, mode="r")
        for read_id in f5file.get_read_ids():
            read = f5file.get_read(read_id)
            bcgrp = read.get_latest_analysis("Basecall_1D") #Basecall_1D_000
            fastq = read.get_analysis_dataset(bcgrp, "BaseCalled_template/Fastq")
            o.out.write(fastq)
        
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    #sys.stderr.write("#Time elapsed: %s\n"%dt)
