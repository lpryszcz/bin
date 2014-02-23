#!/usr/bin/env python
desc="""Annotate DESeq results. 
"""
epilog="""Author: l.p.pryszcz@gmail.com

Mizerow, 2014/02/22
"""

import argparse, gzip, os, sys
from datetime import datetime

def getHandle(handle):
    """Deal with file names and file objects
    Gzipped also supported.
    """
    if type(handle) is str:
        handle = open(handle)
    #open as gzip if needed
    if handle.name.endswith('.gz'):
        handle = gzip.open(handle.name)
    return handle

def load_tab(handle):
    """Return tab-separated file as dict"""    
    data = {}
    for l in getHandle(handle):
        ldata = l[:-1].split('\t')
        if len(ldata) < 2:
            info = "Warning: 2+ fields expected, but %s found! %s\n"
            sys.stderr.write(info % (len(ldata), str(ldata)))
            continue
        #store
        data[ldata[0]] = "\t".join(ldata[1:])
    return data

def deseq2annotate(inputs, annotationFn, pval, verbose):
    """Filter and annotate deseq output"""
    #load annotation
    gene2annotation = load_tab(annotationFn)

    #
    if verbose:
        sys.stderr.write("Processing %s input files...\n" % len(inputs))
    for i, handle in enumerate(inputs, 1):
        #deal with fname instead of fileobj
        handle = getHandle(handle)
        if verbose:
            sys.stderr.write(" %s %s      \r" % (i, handle.name))
        out = open(handle.name+".annotated.txt", "w")
        for l in handle:
            if l.startswith(("#","id\t")):
                out.write(l)
                continue
            ldata = l[:-1].split('\t')
            gene = ldata[0]
            #check pvalue
            p    = ldata[-1]
            if p == "NA":
                continue
            if float(p) >= pval:
                continue
            annotation = ""
            if gene in gene2annotation:
                annotation = gene2annotation[gene]
            out.write("%s\t%s\n"%("\t".join(ldata), annotation))
        out.close()

def main():

    usage  = "%(prog)s [options]"
    parser = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)

    parser.add_argument("-v", "--verbose", default=False, action="store_true")
    parser.add_argument("--version", action="version", version='%(prog)s 0.1')
    parser.add_argument("-i", "--input",      type=file, nargs="+", 
                        help="input streams [stdin]")
    parser.add_argument("-a", "--annotation", type=file, 
                        help="annotation file")
    parser.add_argument("-p", "--pval",       type=float, default=0.05,
                        help="p-value cut-off [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n" % str(o))
    
    deseq2annotate(o.input, o.annotation, o.pval, o.verbose)
        
if __name__=='__main__':
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!                \n")
    dt = datetime.now()-t0
    sys.stderr.write("## Time elapsed: %s\n" % dt)
