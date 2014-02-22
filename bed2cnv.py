#!/usr/bin/env python
desc="""Report CNVs from read counts in genome windows (BED). 
"""
epilog="""Author:
Leszek Pryszcz
l.p.pryszcz@gmail.com

Mizerow, 01/02/2014
"""

import argparse, subprocess, os, sys
import numpy  as np
import pandas as pd
from datetime import datetime

def load_bed(files):
    """Load BED file into Pandas object"""
    for i, file in enumerate(files):
        strain = file.name.split('.')[0] #".".join(file.name.split('.')[:2])
        if i:
            bed = pd.merge(bed, pd.read_table(file, names=('chr', 'start', 'end', strain)), \
                     on=('chr','start','end')) #how='outer'
        else:
            bed = pd.read_table(file, names=('chr','start','end',strain))
                                
    #sort by chromosome position
    bed = bed.sort(columns=('chr','start','end'))
    return bed
    
def bed2cnv(files, out, alpha, verbose):
    """Report deletions/duplications at given threshold"""
    if verbose:
        sys.stderr.write("Loading BED...\n")
    #load BED files
    bed = load_bed(files)

    if verbose:
        sys.stderr.write("Normalising and selecting CNVs...\n")
    #normalise
    deletions = duplications = np.array([False]*len(bed))
    for strain in bed.columns[3:]:
        #reads per 1kb window
        bed[strain] = bed[strain] * 1000.0 / (bed.end - bed.start)
        sys.stderr.write("%s %.2f %.2f\n"%(strain, bed[strain].mean(), bed[strain].std()))
        #observed / expected
        bed[strain] = np.log2(bed[strain] / bed[strain].mean())
        #get deletions
        dels = bed[strain] < bed[strain].quantile(0.0+alpha/2)
        deletions = deletions + dels
        #get duplications
        dups = bed[strain] > bed[strain].quantile(1.0-alpha/2)
        #dups.nonzero()
        duplications = duplications + dups
        sys.stderr.write("%s %s %s\n"%(strain, len(dels.nonzero()[0]), len(dups.nonzero()[0])))

    #select cnvs
    sys.stderr.write("Saving %s deletions and %s duplications.\n" % (deletions.tolist().count(True), duplications.tolist().count(True)))
    bed[duplications + deletions].to_excel(out, sheet_name='CNVs', index=0)
    bed[duplications + deletions].to_csv(out+'.tsv', sep='\t', header=0, index=0)
    
def main():

    usage  = "%(prog)s [options] -i bed1 bed2 bed3"
    parser = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)

    parser.add_argument("-v", "--verbose",  default=False, action="store_true")
    parser.add_argument("--version", action="version", version='%(prog)s 0.1')
    parser.add_argument("-i", "--inputs", nargs="+", type=file,
                        help="input file(s)")
    parser.add_argument("-o", "--output",   default='out.xls', 
                        help="output xls file")
    parser.add_argument("-a", "--alpha",    default=0.05, 
                        help="alpha to call CNVs")

    o = parser.parse_args()
    
    bed2cnv(o.inputs, o.output, o.alpha, o.verbose)

if __name__=='__main__':
    t0=datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!                \n")
    dt=datetime.now()-t0
    sys.stderr.write("## Time elapsed: %s\n" % dt)
