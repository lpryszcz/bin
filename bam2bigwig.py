#!/usr/bin/env python
desc="""Convert BAM to BigWig.
Inspired by: https://www.biostars.org/p/64495/#64680
Added support for non-UCSC genomes.

Require:
- bedGraphToBigWig
- samtools
- pybedtools
- pysam

TBD:
- you can avoid genome fasta if faidx is generated from bam header
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Mizerow, 25/02/2015
"""
import os, pysam, resource, sys
from datetime import datetime
#from pybedtools.contrib.bigwig import bam_to_bigwig
import pybedtools

def get_mapped(bam, verbose=0):
    """Return number of mapped reads in BAM file.
    Create BAM index if not present. """
    # generate BAM index if absent
    if not os.path.isfile(bam+'.bai'):
        cmd = "samtools index %s"%bam
        if verbose:
            sys.stderr.write("[%s]  Indexing BAM file: %s\n"%(datetime.ctime(datetime.now()), cmd))
        os.system(cmd)
        
    # open BAM file
    sam = pysam.AlignmentFile(bam)
    return sam.mapped

def bam2bigwig(bam, genome, output, strand=None, scaled=True, verbose=1):
    """Convert BAM to BigWig scaled in reads per million mapped reads."""
    # skip if outfile exists
    if os.path.isfile(output):
        sys.exit("File exists: %s"%output)
        
    # generate faidx if absent
    faidx = genome+".fai"
    if not os.path.isfile(faidx):
        pysam.FastaFile(genome)
        
    # altered source from https://pythonhosted.org/pybedtools/_modules/pybedtools/contrib/bigwig.html#bam_to_bigwig
    #bam_to_bigwig(bam='path/to/bam', genome='hg19', output='path/to/bigwig')
    if verbose:
        sys.stderr.write("[%s] Converting BAM to BED...\n"%(datetime.ctime(datetime.now()), ))
    kwargs = dict(bg=True, split=True, g=faidx)
    # store strand info
    if strand in ("+", "-", "pos", "neg"):
        if   strand=="pos":
            strand="+"
        elif strand=="neg":
            strand="-"
        kwargs['strand'] = strand
    #store scaling info
    if scaled:
        # speed-up using samtools idxstats
        #readcount = mapped_read_count(bam)
        readcount = get_mapped(bam, verbose)
        _scale = 1 / (readcount / 1e6)
        kwargs['scale'] = _scale
        
    # get genome coverage
    if verbose:
        sys.stderr.write("[%s]  Generating genome coverage\n"%(datetime.ctime(datetime.now()), ))
    x = pybedtools.BedTool(bam).genome_coverage(**kwargs)
    cmds = ['bedGraphToBigWig', x.fn, faidx, output]
    
    # convert to bigWig
    if verbose:
        sys.stderr.write("[%s] Converting BED to bigWig: %s\n"%(datetime.ctime(datetime.now()), " ".join(cmds)))
    os.system(" ".join(cmds))
    
    # clean-up
    os.unlink(x.fn)

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--bam",     required=True,     
                        help="BAM file")
    parser.add_argument("-g", "--genome",  required=True,
                        help="genome FASTA file")
    parser.add_argument("-o", "--output",  required=True,
                        help="output stream   [stdout]")
    parser.add_argument("-s", "--strand", default="both", choices=("both", "+","-", "pos", "neg"), 
                        help="report coverage from + or - strand [%(default)s]")
    parser.add_argument("--scaling",       default=True,  action="store_false",
                        help="disable RPM scaling")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    bam2bigwig(o.bam, o.genome, o.output, o.strand, o.scaling, o.verbose)

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
