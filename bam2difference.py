#!/usr/bin/env python2
desc="""Report sites with differences between bam files. 
"""
epilog="""AUTHOR:
l.p.pryszcz+git@gmail.com

Mizerow, 25/02/2014
"""

import argparse, os, sys
import pysam, subprocess
from datetime import datetime
from bam2snps import _remove_indels

def get_allele_frequency(bases, alphabet="ACGT"):
    """ """
    bases = _remove_indels(bases)
    base2count = {base: 0 for base in alphabet}
    for base in bases.upper():
        if base in base2count:
            base2count[base] += 1

    total_bases = sum(base2count.itervalues())
    if not total_bases:
        return []
    return [(b, c*1.0/total_bases) for b, c in base2count.iteritems()]
            
def bam2difference(fnames, out, minCov, minFreq, homozygous, positions, verbose):
    """
    """
    #write header
    header = "coordinate\t%s\n" % "\t".join(fnames)
    out.write(header)
    
    #open subprocess
    args  = ['samtools', 'mpileup']
    if positions:
        args += ["-l", positions.name]
    args += fnames 
    proc = subprocess.Popen(args, stdout=subprocess.PIPE)

    for line in proc.stdout:
        lData = line.split('\t')
        contig, pos, ref = lData[:3]
        #process all samples
        data = []
        for i in range(3, len(lData), 3):
            cov, bases, quals = lData[i:i+3]
            if int(cov) < minCov:
                break
            #select alleles base on min freq
            alt_alleles = get_allele_frequency(bases)
            alleles = [base for base, freq in filter(lambda x: x[1]>=minFreq, alt_alleles)]
            #skip if not homozygous
            if homozygous and len(alleles)>1:
                break
            #store
            data.append("/".join(sorted(alleles)))
        #report only if more than one allele all samples passed filtering
        if len(set(data))>1 and len(data)==(len(lData)/3)-1:
            out.write("%s:%s\t%s\n" % (contig, pos, "\t".join(data)))
    
def main():
    usage   = "%(prog)s [options]"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
    
    parser.add_argument("-v", "--verbose",  default=False, action="store_true")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument("-i", "--input", nargs="+", 
                        help="input BAM files")
    parser.add_argument("-o", "--output", default=sys.stdout, type=file,
                        help="output stream            [stdout]")
    parser.add_argument("-c", "--cov",   default=10,  type=int,
                        help="min coverage             [%(default)s]")
    parser.add_argument("-f", "--freq",  default=0.2, type=float,
                        help="min allele frequency     [%(default)s]")
    parser.add_argument("--homozygous",  action='store_true', default=False, 
                        help="report only homozygous   [%(default)s]")
    parser.add_argument("--positions",   type=file,
                        help="pre-selected positions   [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write( "Options: %s\n" % str(o) )

    bam2difference(o.input, o.output, o.cov, o.freq, o.homozygous, o.positions, o.verbose)

if __name__=='__main__': 
    t0=datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt=datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
