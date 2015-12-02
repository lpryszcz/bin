#!/usr/bin/env python
desc="""Report genes with evidence of isoform switching. 
For each gene, number of transcripts, cumulative expression,
major isoforms for each condition and each major isoform are reported.
Lowly expressed genes in given condition are coded with `0` (--TPM < 1),
while genes without clear major isoform are marked with `-1` (--frac 0.25).

If you provide transcripts on command-line (-t / --transcripts) the program assumes Salmon output as input.
Otherwise, the cuffdiff output is assumed as input. 


CHANGELOG:
v1.3
- calculate FPKM for salmon output
v1.2
- use TPM for salmon v4.1+ with fpkm column removed
v1.1
- sailfish / salmon supported (.sf)
- support ensembl transcripts fasta files i.e. cdna, ncrna
v1.0
- cufflinks supported (isoforms.fpkm_tracking)
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 6/05/2015
"""

import os, sys
from datetime import datetime
import numpy as np

def get_conditions(line):
    """Return conditions"""
    conditions = []
    lData = line.split('\t')
    for i in range(9, len(lData), 4):
        #strip("_FPKM")
        conditions.append(lData[i][:-5])
    return conditions

def load_gene2transcripts(transcripts):
    """Return gene2transcripts.
    Supports ensembl fasta files.
    """
    # ENSDART00000007748 ensembl:known chromosome:Zv9:21:822304:832471:-1 gene:ENSDARG00000016476
    gene2transcripts = {}
    for line in transcripts:
        if line.startswith(">"):
            lData = line[1:].split()
            tid = lData[0]
            gid = filter(lambda x: x.startswith("gene:"), lData)[0].split(':')[1]
            if gid not in gene2transcripts:
                gene2transcripts[gid] = [tid]
            else:
                gene2transcripts[gid].append(tid)
    return gene2transcripts
        
def parse_fpkm(handle, conditions):
    """Parse fpkm file from cufflinks and yield expression info
    for all trascripts of given gene. """
    pgid = None
    # sort in memory
    data = (l.split('\t') for l in handle)
    for lData in sorted(data, key=lambda x: x[3]): 
        # unload gene/transcript info
        tid, cc, nearest, gid, gene, tss, locus, length, cov = lData[:9]
        # store new gene info
        if pgid != gid:
            if pgid:
                yield pgid, transcripts, fpkms
            # reset
            transcripts = []
            fpkms = [[] for i in range(len(conditions))]
            # store new pgid
            pgid = gid            
        # add transcript info
        transcripts.append(tid)
        # and expression for each condition
        for ci, i in enumerate(range(9, len(lData), 4)):
            fpkm = float(lData[i])
            fpkms[ci].append(fpkm)
    #return conditions, genes, fpkms
    if pgid:
        yield pgid, transcripts, fpkms

def parse_sf(handles, conditions, transcripts):
    """Parse multiple .sf handles and return expression info
    for all trascripts of each gene."""
    # get gene2transcripts
    gene2transcripts = load_gene2transcripts(transcripts)
    # get tid2gid
    tid2gid = {} #tid: gid for
    for gid, tids in gene2transcripts.iteritems():
        for tid in tids:
            tid2gid[tid] = gid
    # prepare
    gene2fpkms = {gid: [[0]*len(tids) for i in range(len(conditions))] \
                  for gid, tids in gene2transcripts.iteritems()}
    # parse handles
    for i, handle in enumerate(handles):
        for line in handle:
            if line.startswith('#'):
                continue
            # unload line
            lData = line[:-1].split('\t')
            if len(lData)==5:
                tid, length, tpm, fpkm, reads = lData
            # salmon v4.1+ https://github.com/COMBINE-lab/salmon/releases
            else:
                tid, length, tpm, reads = lData
            # get gid
            gid = tid2gid[tid]
            # store fpkm for given condition and transcript
            length, tpm = float(length), float(tpm)
            fpkm = tpm*1000.0 / length
            gene2fpkms[gid][i][gene2transcripts[gid].index(tid)] = float(fpkm)
    # yield data
    for gid, fpkms in gene2fpkms.iteritems():
        yield gid, gene2transcripts[gid], fpkms
        
def fpkm2major(handle, out, frac, minTPM, link, transcripts, verbose):
    """Parse expression and report major isoform for each condition"""
    # cufflinks
    if not transcripts:
        # get conditions
        conditions = get_conditions(handle[0].readline())
        if verbose:
            sys.stderr.write("%s conditions: %s\n"%(len(conditions), ", ".join(conditions)))
        parser = parse_fpkm(handle[0], conditions)
    else:
        # get conditions
        ## salmon/RZE024/quant.sf -> RZE024
        conditions = [h.name.split("/")[-2] for h in handle]
        # get parser
        parser = parse_sf(handle, conditions, transcripts)

    # get major isoforms
    header = "#\n#gene id\tno. of transcripts\tFPKM sum\t%s\t%s\n"
    out.write(header%('\t'.join(conditions), '\t'.join(map(str, range(5)))))
    line = "%s\t%s\t%.2f\t%s\t\t%s\n"
    j = k = 0
    for i, (gid, transcripts, TPMs) in enumerate(parser, 1):
        if len(transcripts)<2:
            continue
        k += 1
        #print gid, transcripts, TPMs
        # get max
        t2c = {t: 0 for t in transcripts}
        tmax = []
        for tTPMs in TPMs:
            _tTPMs = sorted(tTPMs, reverse=True)
            # report no major isoform if low expression or small difference
            if   _tTPMs[0] < minTPM:
                major = 0
            elif (1-frac)*_tTPMs[0] < _tTPMs[1]:
                major = -1
            else:
                major = transcripts[tTPMs.index(_tTPMs[0])]
                t2c[major] += 1
            # 
            tmax.append(major)
        # get transcripts that appear most commonly as major
        stranscripts = sorted(t2c.iteritems(), key=lambda x: x[1]>0, reverse=1)
        mtranscripts = [t for t, c in filter(lambda x: x[1], stranscripts)]
        if len(mtranscripts)<2:
            continue
        j += 1
        # recode major isoforms as int
        majors = []
        for m in tmax:
            if m>1:
                majors.append(mtranscripts.index(m)+1)
            else:
                majors.append(m)
        # report
        tpmSum = sum(sum(e) for e in TPMs)
        out.write(line%(_link(link, gid), len(transcripts), tpmSum, "\t".join(map(str, majors)), "\t".join(mtranscripts)))
    if verbose:
        info = "%s genes\n %s genes with 2+ transcripts\n %s genes with 2+ major isoforms\n"
        sys.stderr.write(info%(i, k, j))

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
  
    parser.add_argument('--version', action='version', version='1.3a')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--input", default=[sys.stdin], type=file, nargs="+", 
                        help="isoforms.fpkm_tracking or .sf file(s) from salmon [stdin]")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-f", "--frac", default=0.25, type=float, 
                        help="major isoform has to be larger at least -f than the second most expressed [%(default)s]")
    parser.add_argument("-m", "--minFPKM", default=1.0, type=float, 
                        help="min FPKM to report [%(default)s]")
    parser.add_argument("-t", "--transcripts", type=file, default=None,
                        help="transcripts file in .fasta; needed to get gene2transcripts for .sf input [%(default)s]")
    parser.add_argument("--link", default='=hyperlink("http://www.ensembl.org/Danio_rerio/Gene/Summary?db=core;g=%s", "%s")',
                        help="add hyperlink [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    fpkm2major(o.input, o.output, o.frac, o.minFPKM, o.link, o.transcripts, o.verbose)
 
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror)+str(e)+"\n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
