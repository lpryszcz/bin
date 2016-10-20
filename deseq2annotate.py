#!/usr/bin/env python
desc="""Annotate DESeq results. 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com

Mizerow, 2014/02/22
"""

import argparse, gzip, os, sys, sqlite3
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
   
def load_tab(handles, go):
    """Return tab-separated file as dict"""
    cur = sqlite3.connect(go).cursor()
    data = {}
    for i, handle in enumerate(handles):
        for l in getHandle(handle):
            ldata = l[:-1].split('\t')
            if len(ldata) < 2:
                info = "Warning: 2+ fields expected, but %s found! %s\n"
                sys.stderr.write(info % (len(ldata), str(ldata)))
                continue
            #
            if "GO:" in l:
                for goField in filter(lambda x: "GO:" in x, ldata):
                    cur.execute("select term_type, name from go where acc in ('%s')"%goField.replace(";",",").replace(",","','"))
                    names = "; ".join(name for reltype, name in cur.fetchall())
                    ldata[ldata.index(goField)] = names
            #store
            geneid = ldata[0]
            if geneid not in data:
                data[geneid]  = ["\t".join(ldata[1:])]
            elif len(data[geneid])>i:
                continue
            else:
                data[geneid].append("\t".join(ldata[1:]))
    return data

def to_float(s):
    try:
        return float(s)
    except:
        return s
    
def deseq2annotate(inputs, annotationFns, pvalTh, fcTh, go, verbose):
    """Filter and annotate deseq output"""
    #load annotation
    gene2annotation = load_tab(annotationFns, go)

    #
    if verbose:
        sys.stderr.write("Processing %s input files...\n" % len(inputs))
    gene2data = {}
    samples = []
    for i, handle in enumerate(inputs, 1):
        #deal with fname instead of fileobj
        handle = getHandle(handle)
        samples.append(handle.name)
        if verbose:
            sys.stderr.write(" %s %s      \r" % (i, handle.name))
        for l in handle:
            if l.startswith(("#","id\t")):
                #out.write(l)
                continue
            #unload line
            ldata = l[:-1].split('\t')
            gene = ldata[0]
            baseMean, baseMeanA, baseMeanB, foldChange, log2FoldChange, pval, padj = map(to_float, ldata[1:])
            if gene not in gene2data:
                gene2data[gene] = []
            gene2data[gene].append((baseMeanA, baseMeanB, log2FoldChange, padj))

    #write header
    header=["#gene_id"]
    for sample in samples:
        header.append("control\t%s\tpadj"%sample)
    header.append("annotations\n")
    sys.stdout.write("\t".join(header))
    for gene, edata in gene2data.iteritems():
        #skip if all not significant
        if not filter(lambda x: x[-1]<pvalTh and 0 in x[:2] or \
                      x[-1]<pvalTh and max(x[:2])/min(x[:2])>fcTh, edata):
            continue
        #print edata
        annotation = ""
        if gene in gene2annotation:
            annotation = gene2annotation[gene]
        expression = "\t".join("\t".join(map(str, (baseMeanA, baseMeanB, padj))) \
                               for baseMeanA, baseMeanB, log2FoldChange, padj in edata)
        sys.stdout.write("%s\t%s\t%s\n"%(gene, expression, "\t".join(annotation)))

def main():

    usage  = "%(prog)s [options]"
    parser = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)

    parser.add_argument("-v", "--verbose", default=False, action="store_true")
    parser.add_argument("--version", action="version", version='%(prog)s 0.1')
    parser.add_argument("-i", "--input",      type=file, nargs="+", 
                        help="input streams [stdin]")
    parser.add_argument("-a", "--annotation", nargs="+", type=file, 
                        help="annotation file(s)")
    parser.add_argument("-p", "--pval",       type=float, default=0.05,
                        help="p-value cut-off [%(default)s]")
    parser.add_argument("-f", "--fc",         type=float, default=1.5,
                        help="fold change cut-off [%(default)s]")
    parser.add_argument("--go", default="/home/lpryszcz/cluster/db/gene_ontology.db3",
                        help="GO db3 [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n" % str(o))
    
    deseq2annotate(o.input, o.annotation, o.pval, o.fc, o.go, o.verbose)
        
if __name__=='__main__':
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!                \n")
    dt = datetime.now()-t0
    sys.stderr.write("## Time elapsed: %s\n" % dt)
