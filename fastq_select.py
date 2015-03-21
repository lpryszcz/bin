#!/usr/bin/env python
desc="""Select subset of fastq files.
Support single and paired-end and mate-pair libraries. 
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 6/03/2015
"""

import gzip, os, subprocess, sys
from datetime import datetime
from Bio import SeqIO
from itertools import izip

def get_decompressor(fname, cmd="zcat"):
    """Return subprocess zcat pipe for faster gzip decompression"""
    zcat = subprocess.Popen(['zcat', fname], bufsize=-1, \
                            stdout=subprocess.PIPE)
    # for even more read-intensive tasks, use parallel gzip: pigz
    #zcat = subprocess.Popen(['pigz', '-dc', fname], bufsize=-1, \
    return zcat.stdout

def parse_fastq(handle):
    """Return fastq records as tuples"""
    fastq = []
    for i, l in enumerate(handle, 1):
        fastq.append(l.strip())
        if not i%4:
            yield fastq
            fastq = []

def get_rname(fastq):
    """Return read name from fastq entry"""
    return fastq[0][1:].split()[0]
            
def fastq_select(outbase, files, rnames, verbose):
    """Return reads with given names"""
    # open gzipped files
    if files[0].name.endswith('.gz'):
        files = [get_decompressor(f.name) for f in files]
    # read rnames
    names = set(rname.strip() for rname in rnames if rname.strip())
    if verbose:
        sys.stderr.write(" %s read names loaded...\n"%len(names))
    # init output
    outs = [open("%s.%s.fq"%(outbase, i), "w") for i in range(1, len(files)+1)]
    # process reads
    k = 0
    #parsers = [SeqIO.parse(f, 'fastq') for f in files]
    #for i, mates in enumerate(zip(parsers), 1):
    #for i, mates in enumerate(izip(SeqIO.parse(files[0], 'fastq'), SeqIO.parse(files[1], 'fastq'))):
    for i, mates in enumerate(izip(parse_fastq(files[0]), parse_fastq(files[1]))):
        if verbose and not i%1e5:
            sys.stderr.write(" %i %i \r"%(i, k))
        #if mates[0].name in names:
        if get_rname(mates[0]) in names:
            k += 1
            for r, out in zip(mates, outs):
                #out.write(r.format('fastq'))
                out.write("\n".join(r)+'\n')
            
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--input", nargs="+", type=file, 
                        help="input FASTQ file(s)")
    parser.add_argument("-r", "--rnames", type=file, 
                        help="file with read names to select")
    parser.add_argument("-o", "--outbase", default="out", 
                        help="output base name [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    fastq_select(o.outbase, o.input, o.rnames, o.verbose)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        #sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror)+str(e)+"\n")
        sys.stderr.write("%s\n"%e)
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
