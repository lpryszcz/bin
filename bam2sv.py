#!/usr/bin/env python
desc="""Scan alignments (SAM/BAM) for structural variants.

Detection of deletions, duplications and inversions from paired reads is implemented.
In addition, deletions and duplications are detected from deviations from mean depth
of coverage.

By default, the program dumps reads of interest and depth of coverage information. This
speeds up recalculation by the factor of 20X and should take <1% of BAM file size.

To be done:
+ SV detection
-- insertions testing
-- translocations
-- inversions
+ split read alignment
+ rlen learning
"""
epilog="""Author: l.p.pryszcz@gmail.com
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
        #set variables
        self._set_variables(kwargs)
        #open handles
        self._prepare_handles(bam)
        #prepare storages
        self._init_storages()

    def _set_variables(self, kwargs):
        """Set parameters"""
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
        #set read length
        if 'rlen' in kwargs:
            self.rlen = kwargs['rlen']
        else:
            self.rlen    = None
        #prepare logging
        if   'log' in kwargs:
            self.log = kwargs['log']
        elif 'verbose' in kwargs and kwargs['verbose']:
            self.log = sys.stderr
        else:
            self.log     = None
        #no dump
        if   'nodump' in kwargs:
            self.nodump = kwargs['nodump']
        else:
            self.nodump = False
        #merge by depth of coverage variants
        if   'merge' in kwargs:
            self.merge = kwargs['merge']
        else:
            self.merge = False
        #out
        if 'out' in kwargs:
            self.out = kwargs['out']
        else:
            self.out = sys.stdout
        #coverage fraction
        if 'cov_frac' in kwargs:
            self.cov_frac = kwargs['cov_frac']
        else:
            self.cov_frac = 0.75
        #min dup
        if 'dup_isize_frac' in kwargs:
            self.dup_isize_frac = kwargs['dup_isize_frac']
        else:
            self.dup_isize_frac = 0.9
        #min cnv size from depth of coverage
        if 'cnv_size' in kwargs:
            self.cnv_size = kwargs['cnv_size']
        else:
            self.cnv_size = 1000
        #min cnv size from depth of coverage
        if 'w' in kwargs:
            self.w = kwargs['w']
        else:
            self.w = 100

    def _prepare_handles(self, bam):
        """Open BAM file for reading and set scanning parameters."""
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
        self.del_isize = self.isize_mean + 2*self.isize_stdev
        self.ins_isize = self.isize_mean - 2*self.isize_stdev
        if self.log:
            self.log.write(" FF/FR/RF/RR: %s/%s/%s/%s\n" % tuple(self.pairs))
            self.log.write("  %s chosen\n" % self.orientations[self.orientation])
            info = " median: %.2f  mean: %.2f +- %.2f\n"
            self.log.write(info%(self.isize_median, self.isize_mean, self.isize_stdev))
            info = "  deletion: isize >%.2f\n  insertion: isize <%.2f\n"
            self.log.write(info%(self.del_isize, self.ins_isize))
            
    def _init_storages(self):
        """Initialase storages."""
        #reads storages
        self.delReads = []
        self.dupReads = [] 
        self.insReads = []
        self.invReads = []
        self.traReads = []
        #variants storages
        self.dels = [[] for r in self.refs]
        self.dups = [[] for r in self.refs]
        self.inss = [[] for r in self.refs]
        self.invs = [[] for r in self.refs]
        self.tras = [[] for r in self.refs]

    def alg2orientation(self, alg):
        """Return pair orientation: FF: 0; FR: 1; RF: 2; RR: 4."""
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
        
    def get_isize_stats(self, limit=1e5): 
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
        #update coverage - read count rather that depth of coverage!
        if not os.path.isfile(self.bamdump):
            #self.chr2cov[alg.rname][alg.pos:alg.pos+alg.rlen] += 1
            self.chr2cov[alg.rname][alg.pos] += 1
        #skip right alignment and low quality
        if alg.isize<0 or alg.mapq<self.mapq: # or alg.mrnm < alg.rname
            return
        orient = self.alg2orientation(alg)
        ##insertion/deletion
        #correct pairing
        if  orient == self.orientation:
            ##deletion if significantly larger distance
            if   alg.isize > self.del_isize:
                self.delReads.append(alg)
            ##insertion if significantly smaller distance
            elif alg.isize < self.ins_isize:
                self.insReads.append(alg)
        ##segmental duplication
        #RF <--> FR or FF <--> RR??
        elif self.orientation in (1, 2) and orient in (1, 2) or \
             self.orientation in (0, 4) and orient in (0, 4): 
            self.dupReads.append(alg)
        ##inversion
        #FR/RF -> FF/RR or FF/RR --> FR/RF
        elif self.orientation in (1, 2) and orient in (0, 4) or \
           self.orientation in (0, 4) and orient in (1, 2):
            self.invReads.append(alg)
        ##translocation -- note, some of putative deletions may be also translocations
        #orientation dosn't matter
        if alg.rname != alg.mrnm:
            self.traReads.append(alg)
            
    def sv2bam(self):
        """Dump all alignments important for SVs"""
        if self.log:
            self.log.write("Dumping info to: %s ...\n"%self.bamdump)
        #open out sam/bam handle
        header = self.sam.header
        nalgs = len(self.delReads) + len(self.dupReads) + len(self.insReads) + len(self.invReads) + len(self.traReads)
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
        for algs in (self.delReads, self.dupReads, self.insReads, self.invReads, self.traReads):
            for alg in algs:
                out.write(alg)
        out.close()

    def get_clusters(self, algs):
        """Return clustered algs."""
        #collapse dels by chromosome
        #chr2dels = {i: [] for i, ref in enumerate(self.refs)}
        #py2.6 compatible
        chr2dels = {}
        for i, ref in enumerate(self.refs):
            chr2dels[i] = []
        for alg in algs:
            chr2dels[alg.rname].append(alg)
        clusters = []
        #process each chromosome
        for chri in chr2dels:
            clusters.append([])
            hist = np.zeros((self.lengths[chri]/self.w)+1, dtype=int)
            for alg in chr2dels[chri]:
                hist[alg.pos/self.w] += 1
            #get peaks
            peaks = self.get_peaks(hist, 2*self.w)
            if not peaks:
                continue
            #generate clusters
            i = 0
            for alg in chr2dels[chri]:
                #before current peak
                if alg.pos < peaks[i][0] - self.rlen:
                    continue
                #after current peak
                elif alg.pos > peaks[i][1]:
                    #skip peaks until next is after current read
                    while i < len(peaks) and alg.pos > peaks[i][1]:
                        i += 1
                    if i + 1 >= len(peaks):
                        break
                    #add fresh cluster
                    clusters.append([])
                #store alg to cluster if within peak
                if peaks[i][0] - self.rlen <= alg.pos <= peaks[i][1]:
                    clusters[-1].append(alg)
                    
        #filter by min reads
        clusters = filter(lambda x: len(x) > self.minReads, clusters)#; print len(clusters)
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
        #rlen    = np.mean([alg.rlen for alg in algs])
        #get chromosome info
        chrnames  = [alg.rname for alg in algs]
        mchrnames = [alg.mrnm  for alg in algs]
        return isizes, starts, mstarts, chrnames, mchrnames

    def cnvs_from_pairs(self, reads, storage, cnvType, m=1):
        """Call deletions for paired reads. """
        if not reads:
            return
        #get read clusters
        rlen = self.rlen
        for algs in self.get_clusters(reads):
            isizes, starts, mstarts, chrnames, mchrnames = self.get_algs_features(algs)
            #correct by the read length
            try:
                chri    = int(np.median(chrnames))
            except:
                sys.stderr.write("Cannot get chromosome number: %s\n"%str(chrnames))
                continue
            chrname = self.refs[chri]
            #get left/right mate starts
            leftSt  = np.median(starts)
            rightSt = np.median(mstarts)
            if rightSt<leftSt:
                leftSt, rightSt = rightSt, leftSt
            #check local depth and if enough reads
            nreads  = len(algs)         
            minLcCov = self.w * self.cov_frac * m * self.chr2cov[chri][min(starts):max(starts)+1].mean()
            if nreads < minLcCov:
                continue
            #if too high coverage for deletion
            if   cnvType == "DEL":
                start   = int(leftSt  + self.isize_mean/2.0)
                if start<0: start = 0
                end     = int(rightSt - self.isize_mean/2.0 + rlen) + 1
                size    = end-start
                #check coverage difference - adjust by read start
                cov_obs = self.chr2cov[chri][start-rlen/2:end-rlen/2].mean()
                cov_ratio = cov_obs / self.cov_mean
                if cov_ratio > 1 - self.covD:
                    continue
            #if too low coverage for duplication
            elif cnvType == "DUP":
                start   = int(leftSt  - self.isize_mean/2.0 + rlen)
                if start<0: start = 0
                end     = int(rightSt + self.isize_mean/2.0) + 1
                size    = end-start
                #check dup size
                if size < self.dup_isize_frac*self.isize_mean:
                    continue
                #check coverage difference
                cov_obs = self.chr2cov[chri][start:end].mean()
                cov_ratio = cov_obs / self.cov_mean
                if cov_ratio < 1 + self.covD:
                    continue
            #insertion
            else:
                start = int(leftSt  + self.isize_mean/2.0 - rlen/2)
                if start<0: start = 0
                end   = start + rlen/2
                size  = int(self.isize_mean - np.median(isizes))
                cov_obs = self.chr2cov[chri][start:end].mean()
                cov_ratio = cov_obs / self.cov_mean
                #print min(starts), max(mstarts), leftSt, rightSt, cov_obs, cov_ratio 
                if not 1-self.covD < cov_ratio < 1+self.covD:
                    continue
            #define ploidy
            ploidy  = self.ploidy * cov_ratio
            if end<=start:
                info = "[Warning] End before start: \n %s:%s-%s reads: %s ploidy: %s\n %s\n %s\n %s\n"
                sys.stderr.write(info%(chrname, start, end, nreads, ploidy, \
                                       str(isizes), str(starts), str(mstarts)))
                continue
            #store del
            storage[chri].append((start, end, nreads, ploidy, size))

    def get_peaks(self, hist, offset, th=None):
        """Return collapsed peaks."""
        if th==None:
            th = np.mean(hist) + np.std(hist)
        #collapse neighbours
        peaks = []
        for p, v in enumerate(hist):
            if v<th:
                continue
            #adjust with window
            s = p*self.w
            e = s+self.w
            #extend previous peak if overlap
            if peaks and s <= pe+offset:
                peaks[-1][1] = e
            else:
                peaks.append([s, e])
            pe = e
        return peaks
        
    def _cov2cnv(self, chri, covHist, cnvType, storage, m):
        """Computes CNVs from depth of coverage."""
        th = (m + self.covD*2) * self.cov_mean
        for start, end in self.get_peaks(covHist*m, 2*self.w, th):
            chrname   = self.refs[chri]
            size      = end-start
            nreads    = 0
            if size < self.cnv_size:
                continue
            cov_obs   = self.chr2cov[chri][start:end].mean()
            cov_ratio = cov_obs / self.cov_mean
            ploidy    = self.ploidy * cov_ratio#; print cov_obs, cov_ratio, ploidy
            # check read depth
            if   cnvType=="DUP":
                if cov_ratio < 1+self.covD:
                    continue
            # cnvType=="DEL" and
            elif cov_ratio > 1-self.covD:
                continue
            #skip if already in storage
            overlapping = filter(lambda x: x[0]<start<x[1] or x[0]<end<x[1] \
                                 or start<x[0]<end or start<x[1]<end, \
                                 storage[chri])
            if overlapping:
                txt = "already reported as"
                if self.merge:
                    txt = "replaced"
                    #get positions of overlapping elements
                    idx = [storage[chri].index(o) for o in overlapping]
                    nreads = sum(o[2] for o in overlapping)
                    idx.sort(reverse=True)
                    #remove these events starting from last
                    for i in idx:
                        storage[chri].pop(i)
                    #store event
                    storage[chri].append((start, end, nreads, ploidy, size))
                if self.log:
                    mtmp = "   %s:%s-%s reads:%s ploidy:%.2f"
                    mstr = "\n".join(mtmp%(chrname, o[0], o[1], o[2], o[3]) \
                                        for o in overlapping)
                    info = "  %s %s:%s-%s ploidy:%.2f %s:\n%s\n"
                    self.log.write(info%(cnvType, chrname, start, end, ploidy, txt, \
                                         mstr))
                #skip adding
                continue
            #store del
            storage[chri].append((start, end, nreads, ploidy, size))
            
    def cnvs_from_depth(self):
        """ """
        for chri, chrname in enumerate(self.refs):
            #coverage hist in w size windows
            covHist = np.array([self.chr2cov[chri][i:i+self.w].mean() \
                                for i in range(0, self.lengths[chri]+1, self.w)])
            #duplications and deletions
            for cnvType, storage, m in zip(("DUP", "DEL"), (self.dups, self.dels), \
                                           (1, -1)):
                self._cov2cnv(chri, covHist, cnvType, storage, m)

    def call_variants(self):
        """Call structural variants"""
        #get expected coverage
        self.cov_mean = 1.0 * sum(sum(c) for c in self.chr2cov) / sum(self.lengths)
        self.minReads = self.cov_frac * self.cov_mean * self.w
        if self.log:
            info = "Calling variants...\n Expected read coverage: %.3f\n"
            self.log.write(info % self.cov_mean)
        #write header
        self.out.write("#chrom\tstart\tend\tname\treads pairs\tploidy\tsize\n")
        #set output formats
        self.cnvline = "%s\t%s\t%s\t%s\t%s\t%.2f\t%s\n"
        #call deletions and duplications
        for reads, storage, cnvType in zip((self.delReads, self.dupReads), \
                                           (self.dels, self.dups), ("DEL", "DUP")):
            #call from paired reads
            self.cnvs_from_pairs(reads, storage, cnvType)
        #call from read depth
        self.cnvs_from_depth()
        #insertions
        self.cnvs_from_pairs(self.insReads, self.inss, "INS", 4)
        
        ##report
        for sv, cnvType in zip((self.dels, self.dups, self.invs), \
                               ("DEL", "DUP", "INS")):
            i = 0
            for chri, variants in enumerate(sv):
                variants.sort()
                for i, (start, end, nreads, ploidy, size) in enumerate(variants, i+1):
                    chrname = self.refs[chri]
                    #define name
                    name    = "%s%3i" % (cnvType, i)
                    name    = name.replace(" ","0")
                    #write output
                    self.out.write(self.cnvline%(chrname, start, end, name, nreads, \
                                                 ploidy, size))
            if self.log:
                sys.stderr.write(" %ss: %s\n" % (cnvType, i))
                
                
        
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
                self.log.write(info % (i, i*100.0/self.nalgs, '%', len(self.delReads), \
                                len(self.dupReads), len(self.insReads), len(self.invReads), len(self.traReads), \
                                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
            #skip unmapped and secondary alignments
            if alg.rname<0 or alg.is_secondary:
                continue
            #add read
            self.add_read(alg)
        if self.log:
            self.log.write(" %s alignments parsed. \n"%i)
        #dump all important info
        if not self.nodump and not os.path.isfile(self.bamdump):
            self.sv2bam()
        #get mean rlen
        if not self.rlen:
            self.rlen = np.mean([alg.rlen for alg in self.delReads])
            if self.log:
                self.log.write(" Mean read length: %.2f \n"%self.rlen)
        #call variants
        self.call_variants()        
            
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--bam", required=True,       
                        help="BAM file")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-p", "--ploidy",    default=2, type=int, 
                        help="ploidy          [%(default)s]")
    parser.add_argument("-q", "--mapq",      default=20, type=int, 
                        help="min mapping quality for variants [%(default)s]")
    parser.add_argument("--rlen",            default=None, type=int, 
                        help="read length     [get from data]")
    parser.add_argument("-c", "--covD",      default=0.33, type=float, 
                        help="min coverage change to call deletion/duplication [%(default)s]")
    parser.add_argument("--cov_frac",        default=0.1, type=float, 
                        help="min fraction of local depth to call variation [%(default)s]")
    parser.add_argument("--dup_isize_frac",  default=0.9, type=float, 
                        help="min duplication size as insert size fraction [%(default)s]")
    parser.add_argument("--cnv_size",        default=1000, type=int, 
                        help="min CNV size from depth of coverage [%(default)s]")
    parser.add_argument("--merge",           default=False,  action="store_true",
                        help="merge read pairs variants using depth of coverage variants")
    parser.add_argument("--nodump",          default=False,  action="store_true",
                        help="dump SV reads for faster recalculations")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    #initialise structural variants
    sv = SVs(o.bam, out=o.output, mapq=o.mapq, ploidy=o.ploidy, covD=o.covD, \
             cov_frac=o.cov_frac, rlen=o.rlen, dup_isize_frac=o.dup_isize_frac, \
             cnv_size=o.cnv_size, merge=o.merge, \
             nodump=o.nodump, verbose=o.verbose)
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
