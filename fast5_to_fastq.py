#!/usr/bin/env python3
desc="""Report FastQ from basecalled Fast5 file(s).

Dependencies: h5py

TO DO:
- 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Barcelona, 13/08/2019
"""

import h5py, os, sys
from datetime import datetime

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='0.10a')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument("-i", "--fast5", nargs="+", help="input Fast5 file(s)")
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("w"), help="output stream [stdout]")
    parser.add_argument("--basecall_group", default="Basecall_1D_000", help="basecall group to use in Fast5 file")
    parser.add_argument("--replace", action='store_true', help="process also reads with Analyses/Resquiggle_* entries")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    for fn in o.fast5:
        h5 = h5py.File(fn, 'r')
        for i, rname in enumerate(h5):
            #if i>1: break
            # skip entries that are already processed
            resquiggle_group = "Resquiggle_1D_%s"%o.basecall_group.split("_")[-1]
            if resquiggle_group in h5["%s/Analyses"%rname].keys() and not o.replace:
                break
            fastq = h5["%s/Analyses/%s/BaseCalled_template/Fastq"%(rname, o.basecall_group)][()].tostring().decode()
            o.out.write(fastq)
        h5.close()
        
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
