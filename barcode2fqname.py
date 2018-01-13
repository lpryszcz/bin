#!/usr/bin/env python
# Process PE reads and insert barcode of given size (size) into fq name. 
# New files are written in current directory! 
# USAGE: barcode2fqname.py size reads_?.fq.gz

import os, sys, gzip, subprocess
from itertools import izip as zip

def fastq_iterator(handle):
    """Return FastQ fields"""
    data = []
    for i, l in enumerate(handle):
        if not i%4:
            if len(data)==4: yield data
            data = []
        data.append(l[:-1])
    if len(data)==4: yield data

def barcode2fqname(size, fnames):
    """Process barcoded FastQ and insert barcode into name field"""
    barcodes = {}
    i = 0
    handle1, handle2 = [subprocess.Popen(['zcat', fn], stdout=subprocess.PIPE).stdout for fn in fnames]
    #handle1, handle2 = [gzip.open(fn) for fn in fnames]
    out1, out2 = [gzip.open(fn, "w") for fn in map(os.path.basename, fnames)]
    for i, (fq1, fq2) in enumerate(zip(fastq_iterator(handle1), fastq_iterator(handle2))):
        if not i%1e3: sys.stderr.write(' %s %s\r'%(i, len(barcodes)))
        barcode = fq1[1][:size]
        if barcode not in barcodes:
            barcodes[barcode] = 1
        else:
            barcodes[barcode] += 1
        out1.write("@%s.%s/1\n%s\n+\n%s\n"%(barcode, barcodes[barcode], fq1[1][size:], fq1[3][size:]))
        out1.write("@%s.%s/2\n%s\n+\n%s\n"%(barcode, barcodes[barcode], fq2[1], fq2[3]))
    sys.stdout.write("%s barcodes in %s sequences\n"%(len(barcodes), i))

if __name__=="__main__": 
    barcode2fqname(int(sys.argv[1]), sys.argv[2:4])
