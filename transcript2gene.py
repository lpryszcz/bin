#!/usr/bin/env python2
desc="""Combine values for transcripts and report summed values for genes.

CHANGELOG:
- v1.2
-- report also gene_name beside geneid
-- added hyperlinks
- v1.1
-- .tsv annotation support
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Warsaw, 30/09/2015
"""

import os, sys
import numpy as np
from datetime import datetime

def get_transcript2gene_tsv(handle):
    """Load transcript-gene relationships from tab-delimited file."""
    tid2gid = {}
    for l in handle:
        tid, gid = l[:-1].split()
        tid2gid[tid] = gid
    return tid2gid
    
def get_transcript2gene_gtf(handle):
    """Load transcript-gene relationships from ensembl GTF."""
    tid2gid = {}
    for l in handle:
        l = l.strip()
        if l.startswith('#') or not l: 
            continue
        
        contig,source,feature,start,end,score,strand,frame,comments = l.split('\t')
        if feature != "transcript":
            continue
    
        description={}
        for atr_value in comments.split(';'):
            atr_value = atr_value.strip()
            if not atr_value:
                continue
            atr   = atr_value.split()[0]
            value = " ".join( atr_value.split()[1:] ).strip('"')
            #value = value.strip('"')
            description[atr]=value

        #
        if "transcript_id" in description and "gene_id" in description:
            tid = description["transcript_id"]
            gid = description["gene_id"]
            genename = "-"
            if "gene_name" in description:
                genename = description["gene_name"]
            tid2gid[tid] = (gid, genename)
    return tid2gid
   
def transcript2gene(handle, out, gtf, tsv, header=0, link="", verbose=0):
    """Report summed gene expression from transcripts"""
    if verbose:
        sys.stderr.write("Parsing annotation...\n")
    if gtf:
        tid2gid = get_transcript2gene_gtf(gtf)
    elif tsv:
        tid2gid = get_transcript2gene_tsv(tsv)
        
    if verbose:
        sys.stderr.write(" %s transcripts parsed.\n"%len(tid2gid))

    if verbose:
        sys.stderr.write("Parsing input...\n")
    gid2data = {}
    for i, l in enumerate(handle):
        # write header
        lData = l[:-1].split('\t')
        if i<header:
            if gtf:
                lData = ["geneid", "genename"] + lData[1:]
            elif tsv:
                lData = ["geneid", ] + lData[1:]
            out.write("\t".join(lData)+"\n")
            continue
        # unload data
        tid, info = lData[0], map(float, lData[1:])
        # check if tid in dict
        if tid not in tid2gid:
            sys.stderr.write("[WARNING] Transcript '%s' not found in tid2gid!\n"%tid)
            continue
        gid = tid2gid[tid]
        # add gid
        if gid not in gid2data:
            gid2data[gid] = []
        # store info
        gid2data[gid].append(info)
    # sum up
    for gid, values in gid2data.iteritems():
        geneid = gid
        if gtf:
            geneid, genename = gid
        # sum values for each column
        a = np.array(values)
        summed = "\t".join(map(str, a.sum(axis=0)))
        if gtf:
            out.write("%s\t%s|%s\t%s\n"%(_link(link, geneid), geneid, genename, summed))
        elif tsv:
            out.write("%s\t%s\n"%(_link(link, geneid), summed))

def _link(link, gid):
    """Return link"""
    if link:
        return link % tuple([gid]*link.count('%s'))
    return gid
        
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.2b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--input", default=sys.stdin, type=file,  
                        help="input stream    [stdin]")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("--header", default=0, type=int, 
                        help="header lines [%(default)s]")
    parser.add_argument("--link", default='=hyperlink("http://www.ensembl.org/Danio_rerio/Gene/Summary?db=core;g=%s", "%s")',
                        help="add hyperlink [%(default)s]")                        
    # mutually exclusive annotation
    annota = parser.add_mutually_exclusive_group(required=True)
    annota.add_argument("-g", "--gtf",   type=file, 
                        help="annotation .gtf")
    annota.add_argument("-t", "--tsv",   type=file, 
                        help="annotation .tsv")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    transcript2gene(o.input, o.output, o.gtf, o.tsv, o.header, o.link, o.verbose)

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
