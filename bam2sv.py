#!/usr/bin/env python
desc="""Scan alignments (BAM) for structural variants.

TBD:
- deletions: based on depth of coverage and reads pairing
-- micro-deletions: instability sites
- duplications
- translocations
- inversions
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow, 13/03/2014
"""

import os, sys
import pickle, pysam, resource
from datetime import datetime
import numpy as np
from scipy import stats, signal

class SVs(object):
    """Store BAM related data. Take BAM file name as input."""
    def __init__(self, bam, **kwargs): 
        #define mapq
        if 'mapq' in kwargs:
            self.mapq = kwargs['mapq']
        else:
            self.mapq = 20
        #ploidy
        if 'ploidy' in kwargs:
            self.ploidy = kwargs['ploidy']
        else:
            self.ploidy = 2
        #q - percentile
        if 'q' in kwargs:
            self.q = kwargs['q']
        else:
            self.q = 1.0
        #min coverage change
        if 'covD' in kwargs:
            self.covD = kwargs['covD']
        else:
            self.covD = 0.33
        #min coverage change
        if 'rlen' in kwargs:
            self.rlen = kwargs['rlen']
        else:
            self.rlen    = 100
        #prepare logging
        if   'log' in kwargs:
            self.log = kwargs['log']
        elif 'verbose' in kwargs and kwargs['verbose']:
            self.log = sys.stderr
        else:
            self.log     = None
        #out
        if 'out' in kwargs:
            self.out = kwargs['out']
        else:
            self.out = sys.stdout
        #define internal sam and shortcuts
        self.bam     = bam
        self.bamdump = self.bam + ".sv_algs.bam"
        if os.path.isfile(self.bamdump):
            if self.log:
                self.log.write("Loading data from dump: %s ...\n" % self.bamdump)
            self.sam     = pysam.Samfile(self.bamdump)
            self.chr2cov = pickle.load(open(self.bamdump+'.chr2cov.pickle'))
            self.isize_mean   = float(self.sam.header['PG'][0]['VN'])
            self.isize_stdev  = float(self.sam.header['PG'][1]['VN'])
            self.isize_median = float(self.sam.header['PG'][2]['VN'])
            self.pairs        = [int(x) for x in self.sam.header['PG'][3]['VN'].split(",")]
            self.nalgs        = int(self.sam.header['PG'][4]['VN'])
        else:
            self.sam   = pysam.Samfile(self.bam)
            self.nalgs = self.sam.mapped + self.sam.unmapped
        #prepare storages
        self.reset()
        self.ndels = self.nins = self.ndups = self.ninvs = self.ntrans = 0
        #shortcuts
        self.refs    = self.sam.references
        self.lengths = self.sam.lengths
        #coverage data storage
        if not os.path.isfile(self.bamdump):
            self.chr2cov = [np.zeros(l, dtype='int') for l in self.lengths]
            #estimate insert size statistics
            self.get_isize_stats()
        #select main orientation
        self.orientations = ("FF", "FR", "RF", "RR")
        self.orientation  = self.pairs.index(max(self.pairs))
        #define deletion and insertion thresholds
        self.del_isize = self.isize_mean + 2* self.isize_stdev
        self.ins_isize = self.isize_mean - 2* self.isize_stdev
        if self.log:
            self.log.write(" FF/FR/RF/RR: %s/%s/%s/%s\n" % tuple(self.pairs))
            self.log.write("  %s chosen\n" % self.orientations[self.orientation])
            self.log.write(" median: %.2f  mean: %.2f +- %.2f\n" % \
                           (self.isize_median, self.isize_mean, self.isize_stdev))
            self.log.write("  deletion: isize >%.2f\n  insertion: isize <%.2f\n" % \
                           (self.del_isize, self.ins_isize))

    def reset(self, isize=False):
        """Reset storage."""
        self.dels    = []
        self.dups    = [] 
        self.ins     = []
        self.invs    = []
        self.trans   = []

    def alg2orientation(self, alg):
        """Return pair orientation.
        FF: 0
        FR: 1
        RF: 2
        RR: 4
        """
        ##FR/RF
        if alg.is_reverse != alg.mate_is_reverse:
            #FR
            if alg.is_read1 and not alg.is_reverse or \
               alg.is_read2 and not alg.is_reverse:
                return 1
            #RF
            else:
                return 2
        #RR - double check that!
        elif alg.is_read1 and alg.is_reverse or \
             alg.is_read2 and not alg.is_reverse:
            return 3
        #FF
        else:
            return 0
        
    def get_isize_stats(self, limit=1e6): 
        """Estimate insert size median, mean and stdev.
        Also count pair orientations and select main.
        """
        if self.log:
            self.log.write("Estimating insert size stats...\n")
        isizes = []
        self.pairs = [0, 0, 0, 0]
        for alg in pysam.Samfile(self.bam):
            #take only reads with good alg quality and one read per pair
            if alg.mapq < self.mapq or alg.isize < 1:
                continue
            #store isize
            isizes.append(alg.isize)
            #store pair orientation
            self.pairs[self.alg2orientation(alg)] += 1
            #stop if limit reached
            if len(isizes) >= limit:
                break
        #get rid of right 5 percentile
        maxins = stats.scoreatpercentile(isizes, 100-self.q)
        minins = stats.scoreatpercentile(isizes, self.q)
        isizes = filter(lambda x: minins<x<maxins, isizes)
        #store
        self.isize_median = np.median(isizes)
        self.isize_mean   = np.mean(isizes)
        self.isize_stdev  = np.std(isizes)
        
    def add_read(self, alg):
        """Update handles for coverage, insert size, etc."""
        #update coverage
        if not os.path.isfile(self.bamdump):
            self.chr2cov[alg.rname][alg.pos:alg.pos+alg.rlen] += 1
        orient = self.alg2orientation(alg)
        ##insertion/deletion
        #correct pairing
        if  orient == self.orientation:
            ##deletion if significantly larger distance
            if   alg.isize > self.del_isize:
                self.dels.append(alg)
            ##insertion if significantly smaller distance
            elif alg.isize < self.ins_isize:
                self.ins.append(alg)
        ##segmental duplication
        #RF <--> FR or FF <--> RR??
        elif self.orientation in (1, 2) and orient in (1, 2) or \
             self.orientation in (0, 4) and orient in (0, 4): 
            self.dups.append(alg)
        ##inversion
        #FR/RF -> FF/RR or FF/RR --> FR/RF
        elif self.orientation in (1, 2) and orient in (0, 4) or \
           self.orientation in (0, 4) and orient in (1, 2):
            self.invs.append(alg)
        ##translocation -- note, some of putative deletions may be also translocations
        #orientation dosn't matter
        if alg.rname != alg.mrnm:
            self.trans.append(alg)
            
    def dump(self):
        """Dump all alignments important for SVs"""
        if self.log:
            self.log.write("Dumping info to: %s ...\n"%self.bamdump)
        #open out sam/bam handle
        header = self.sam.header
        nalgs = len(self.dels) + len(self.dups) + len(self.ins) + len(self.invs) + len(self.trans)
        isize_info = [{"ID": "isize_mean",   "VN": self.isize_mean},
                      {"ID": "isize_stdev",  "VN": self.isize_stdev},
                      {"ID": "isize_median", "VN": self.isize_median},
                      {"ID": "pairs", "PN": "/".join(self.orientations),
                       "VN": ",".join(str(x) for x in self.pairs)},
                      {"ID": "nalgs",  "VN": nalgs}]
        header['PG'] = isize_info 
        out = pysam.Samfile(self.bamdump, "wb", header=header)
        #store per base depth of coverage
        with open(self.bamdump+".chr2cov.pickle", "w") as f:
            pickle.dump(self.chr2cov, f, 2)
        #dump individual algs
        for algs in (self.dels, self.dups, self.ins, self.invs, self.trans):
            for alg in algs:
                out.write(alg)
        out.close()

    def get_clusters(self, algs, w=100):
        """Return clustered algs."""
        #collapse dels by chromosome
        chr2dels = {} #i: [] for i, ref in enumerate(self.refs)}
        for alg in algs:
            if alg.rname not in chr2dels:
                chr2dels[alg.rname] = []
            chr2dels[alg.rname].append(alg)
        clusters = []
        #process each chromosome
        for chromi in chr2dels:
            clusters.append([])
            hist = np.zeros(self.lengths[chromi]/w, dtype=int)
            for alg in chr2dels[chromi]:
                hist[alg.pos/w] += 1
            #get peaks
            peaks = signal.find_peaks_cwt(hist, np.arange(1+(50/w), 200/w), \
                                          min_snr=1, noise_perc=10)
            if not peaks:
                continue
            #TO ADD collapse neighbours
            
            #adjust with window
            peaks = [p*w for p in peaks]
            #generate clusters
            i = 0
            for alg in chr2dels[chromi]:
                #before current peak
                if alg.pos < peaks[i] - w/2 - self.rlen:
                    continue
                #after current peak
                elif alg.pos > peaks[i] + w/2:
                    #skip peaks until next is after current read
                    while i < len(peaks) and alg.pos > peaks[i] + w/2:
                        i += 1
                    if i + 1 >= len(peaks):
                        break
                    #add fresh cluster
                    clusters.append([])
                #store alg to cluster if within peak
                if peaks[i] - w/2 - self.rlen <= alg.pos <= peaks[i] + w/2:
                    clusters[-1].append(alg)
                    
        #filter by min reads
        clusters = filter(lambda x: len(x) > self.minReads, clusters)
        return clusters
        
    def get_algs_features(self, algs):
        """Return algs starts, mate starts, isizes, r"""
        #filter by isize percentile    
        isizes  = [alg.isize for alg in algs]
        min_isize = stats.scoreatpercentile(isizes, self.q)
        max_isize = stats.scoreatpercentile(isizes, 100-self.q)
        algs    = filter(lambda x: min_isize <= x.isize <= max_isize, algs)
        #get sizes etc
        isizes  = [alg.isize for alg in algs]
        starts  = [alg.pos   for alg in algs]
        mstarts = [alg.mpos  for alg in algs]
        rlen    = np.mean([alg.rlen for alg in algs])
        #get chromosome info
        chrnames  = [alg.rname for alg in algs]
        mchrnames = [alg.mrnm  for alg in algs]
        return isizes, starts, mstarts, rlen, chrnames, mchrnames

    def call_deletions(self):
        """Call deletions for paired reads.
        ADD calling from depth-of-coverage. 
        """
        if not self.dels:
            return
        #get read clusters
        clusters = self.get_clusters(self.dels)
        #write header
        self.out.write("##BED-like format\n#chrom\tstart\tend\tname\tscore\tploidy\n")
        bedline = "%s\t%s\t%s\t%s\t%s\t%.2f\n"
        for algs in clusters:
            isizes, starts, mstarts, rlen, chrnames, mchrnames = self.get_algs_features(algs)
            #correct by the read length
            chri    = int(np.median(chrnames))
            chrname = self.refs[chri]
            start   = int(np.mean(starts)  + self.isize_mean/2.0 + rlen)
            end     = int(np.mean(mstarts) - self.isize_mean/2.0 + 2*rlen) + 1
            size    = int(np.mean(isizes)  - self.isize_mean)
            #check coverage difference
            cov_obs = 1.0 * sum(self.chr2cov[chri][start:end]) / size
            cov_ratio = cov_obs / self.cov_mean
            #apparent deletion
            if cov_ratio > 1-self.covD:
                continue
            self.ndels += 1
            name    = "DEL%3i" % self.ndels
            ploidy  = self.ploidy * cov_ratio
            name    = name.replace(" ","0")
            score   = len(algs)
            self.out.write(bedline % (chrname, start, end, name, score, ploidy))
            
    def call_variants(self):
        """Call structural variants"""
        #get expected coverage
        self.cov_mean = 1.0 * sum(sum(c) for c in self.chr2cov) / sum(self.lengths)
        #set variables base on expected coverage
        self.minReads = int(0.1 * self.cov_mean)
        if self.log:
            info = "Calling variants...\n Expected coverage: %.2f\n  >%s reads to call variant\n"
            self.log.write(info%(self.cov_mean, self.minReads))         
        #call deletions
        self.call_deletions()
        #call deletions using depth of coverage approach
        
    def parse(self, test=0):
        """Parse sam alignments and store info"""
        #parse algs
        if self.log:
            self.log.write("Parsing alignments...\n")
        pchrom = ""
        for i, alg in enumerate(self.sam, 1):
            if test and i > test:
                break
            #write log
            if self.log and not i % 1e5:
                info = " %s [%.1f%s]  reads for dels: %s dups: %s ins: %s invs: %s trans: %s [%s Mb]\r"
                self.log.write(info % (i, i*100.0/self.nalgs, '%', len(self.dels), \
                                len(self.dups), len(self.ins), len(self.invs), len(self.trans), \
                                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
            #skip low quality alignments - alg.rname=-1 for unaligned
            #add each pair only once
            #ie transclocations read ref < mate ref
            if alg.isize<0 or alg.mapq<self.mapq or alg.rname<0 or alg.mrnm < alg.rname \
               or alg.is_secondary:
                continue
            #add read
            self.add_read(alg)
        if self.log:
            self.log.write(" %s alignments parsed. \n"%i)
        #dump all important info
        #if not os.path.isfile(self.bamdump):
        #    self.dump()
        #call variants
        self.call_variants()        
            
def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-i", "--bam",       
                        help="BAM file")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-p", "--ploidy",    default=2, type=int, 
                        help="ploidy          [%(default)s]")
    parser.add_argument("-q", "--mapq",      default=10, type=int, 
                        help="min mapping quality [%(default)s]")
    parser.add_argument("--rlen",            default=100, type=int, 
                        help="read length     [%(default)s]")
    parser.add_argument("-c", "--covD",      default=0.33, type=float, 
                        help="min coverage change [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    #initialise structural variants
        sv = SVs(o.bam, out=o.output, mapq=o.mapq, ploidy=o.ploidy, covD=o.covD, \
                 rlen=o.rlen, verbose=o.verbose)
    #call variants in all chromosomes
    sv.parse()

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
