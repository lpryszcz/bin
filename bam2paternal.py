#!/usr/bin/env python
desc="""Calculate paternal expression.

CHANGELOG:
1.0b
- recognise the introns properly
1.0a
- 
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw, 17/11/2015
"""

import argparse, os, re, sys
from datetime import datetime
import subprocess
import numpy as np

read_start_pat = re.compile('\^.')
indel_pat = re.compile('[+-]\d+')

def _remove_indels(alts):
    """
    Remove indels from mpileup.     .$....,,,,....,.,,..,,.,.,,,,,,,....,.,...,.,.,....,,,........,.A.,...,,......^0.^+.^$.^0.^8.^F.^].^],
    ........,.-25ATCTGGTGGTTGGGATGTTGCCGCT..
    And ignore introns:     >>><<<<><><
    
    """
    #But first strip start/end read info.
    alts = "".join(read_start_pat.split(alts)).replace('$', '')
    #remove indels info
    m = indel_pat.search(alts)
    while m:
        #remove indel
        pos = m.end() + int(m.group()[1:])
        alts = alts[:m.start()] + alts[pos:]
        #get next match
        m = indel_pat.search(alts, m.start())
    return alts

def get_alt_allele(base_ref, cov, alg, minFreq, alphabet, reference, bothStrands):
    """Return alternative allele only if different than ref and freq >= minFreq."""
    #remove deletions
    alts = alg
    dels = alts.count('*') 
    #remove insertions
    alts = _remove_indels(alts)
    #get base counts
    baseCounts = [(alts.upper().count(base), base) for base in alphabet]
    #get base frequencies
    for base_count, base in sorted(baseCounts):
        freq = base_count*1.0/len(alts)#cov
        if base!=base_ref and freq >= minFreq:
            #check if alt base in both strands
            if bothStrands: 
                if not base.upper() in alts or not base.lower() in alts:
                    return
            return (base, freq) # base!=base_ref and

def get_major_alleles(cov, alg, minFreq, alphabet, bothStrands):
    """Return major alleles that passed filtering and their frequencies."""
    #remove deletions
    alts = alg
    dels = alts.count('*') 
    #remove insertions
    alts = _remove_indels( alts )
    #get base frequencies
    bases, freqs = [], []
    for base_count, base in zip((alts.upper().count(b) for b in alphabet), alphabet):
        freq = base_count*1.0/cov
        #check freq
        if freq < minFreq:
            continue
        #check if alt base in both strands
        if bothStrands: 
            if not base.upper() in alts or not base.lower() in alts:
                continue
        bases.append(base)
        freqs.append(freq)
    return bases, freqs

def concatenate_mpileup(data):
    """Concatenate mpileup from multiple samples / replicas"""
    cov, algs, quals = 0, "", ""
    for _cov, _algs, _quals in zip(data[0::3], data[1::3], data[2::3]):
        cov += int(_cov)
        if int(_cov):
            algs  += _algs
            quals += _quals
    return cov, algs, quals

def get_meanQ(quals, offset=33):
    """Return mean quality"""
    return np.mean([ord(q)-offset for q in quals])
    
def get_frequency(refBase, bases, freqs):
    """Return frequency of refBase"""
    freq = 0
    if refBase in bases:
        idx = bases.index(refBase)
        freq = freqs[idx]
    return freq
    
def parse_mpileup(bam, minDepth, minfreq, mpileup_opts, regions,
                  report_missing, verbose, bothStrands=0, alphabet='ACGT', missing="-"):
    """Run mpileup subprocess and parse output."""
    # add regions to mpileup opts
    if regions:
        mpileup_opts += ' -l %s' % regions
    #open subprocess
    args = ['samtools', 'mpileup'] + mpileup_opts.split() + bam
    if verbose:
        sys.stderr.write("Running samtools mpileup...\n %s\n" % " ".join(args))
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, bufsize=65536)
    for line in proc.stdout:
        line      = line.strip()
        lineTuple = line.split('\t')
        #get coordinate
        contig, pos, baseRef = lineTuple[:3]
        # load DNAseq / EXOMEseq for female & male, and RNAseq
        fData = lineTuple[3:6]
        mData = lineTuple[6:9]
        rData = lineTuple[9:]

        fCov, fAlgs, fQuals = fData # concatenate_mpileup(fData)
        mCov, mAlgs, mQuals = mData # concatenate_mpileup(mData)
        fCov, mCov = int(fCov), int(mCov)
        if fCov < minDepth or mCov < minDepth:
            continue
        fRef, fFreq = get_major_alleles(fCov, fAlgs, minfreq, alphabet, bothStrands)
        mRef, mFreq = get_major_alleles(mCov, mAlgs, minfreq, alphabet, bothStrands)
        if not fRef or not mRef or set(fRef).intersection(mRef):
            continue
        #print fRef, mRef, fRef==mRef
        # get main base of female and male
        fBase = fRef[0]
        mBase = mRef[0]            
        # mean qual
        fmeanQ = get_meanQ(fQuals)
        mmeanQ = get_meanQ(mQuals)
        # report only freq of maternal and paternal alleles from RNAseq
        sFreqs = []
        mCount = 0
        for cov, alg, quals in zip(rData[0::3], rData[1::3], rData[2::3]):
            cov = int(cov)
            meanQ = get_meanQ(quals)
            if cov<minDepth:
                if report_missing:
                    mCount += 1
                    sFreqs.append("%i\t%s\t%s\t%s"%(cov, meanQ, missing, missing))
                continue
            # check for SNP
            bases, freqs = get_major_alleles(cov, alg, 0.000001, alphabet, bothStrands)
            if not bases or bases==baseRef:
                if report_missing:
                    mCount += 1
                    sFreqs.append("%i\t%s\t%s\t%s"%(cov, meanQ, missing, missing))
                continue
            # get freq of female / male base
            _fFreq = get_frequency(fBase, bases, freqs)
            _mFreq = get_frequency(mBase, bases, freqs)            
            sFreqs.append("%i\t%s\t%s\t%s"%(cov, meanQ, _fFreq, _mFreq))
        # skip if not all samples passsed or all samples missing if report_missing=True
        if len(sFreqs)<len(bam)-2 or mCount==len(bam)-2:
            continue
        yield (contig, pos, fBase, mBase, fCov, fmeanQ, mCov, mmeanQ, '\t'.join(sFreqs))
        
def main():

    usage  = "%(prog)s [options] -b femaleDNA.bam maleDNA.bam rnaSeq*.bam" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0b')
    parser.add_argument("-b", "--bam", nargs="+", 
                        help="input RNA-Seq BAM file(s)")
    parser.add_argument("-l", "--regions", default="", 
                        help="list of positions (chr pos) or regions (BED)")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-r", "--report_missing", default=False, action="store_true",  
                        help="report positions with missing values")
    parser.add_argument("--minDepth", default=5,  type=int,
                        help="minimal depth of coverage [%(default)s]")
    parser.add_argument("--minfreq",  default=0.01, type=float,
                        help="min frequency for base calling [%(default)s]")
    parser.add_argument("--mpileup_opts",   default="-I -q 15 -Q 20",  
                        help="options passed to mpileup         [%(default)s]")
                        
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    if len(o.bam)<3:
        sys.exit("Provide at least 3 BAM files!")
        
    #parse pileup in single-thread mode 
    parser = parse_mpileup(o.bam, o.minDepth, o.minfreq, o.mpileup_opts, o.regions, \
                           o.report_missing, o.verbose)
    # header
    header = "#chrom\tposition\tfemale base\tmale base\tfemale coverage\tfemale meanQ\tmale coverage\tmale meanQ"
    for fn in o.bam[2:]:
        header += "\t%s coverage\t%s meanQ\t%s female base freq\t%s male base freq"%tuple([fn,]*4)
    o.output.write(header+"\n")
    # report 
    info = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
    for data in parser:
        o.output.write(info%data)
    
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
