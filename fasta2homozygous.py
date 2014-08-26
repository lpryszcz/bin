#!/usr/bin/env python
desc="""Align genome onto itself (BLAT) and remove heterozygous (redundant) scaffolds.

TO ADD:
- scaffold extension based on overlapping matches
- reporting of haplotypes
- recognise heterozygous contigs with translocations
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 26/08/2014
"""

import gzip, math, os, sys
from datetime import datetime
from Bio import SeqIO

def blat(fasta, identity, verbose):
    """Start BLAT"""
    #prepare BLAT command
    identity = int(100*identity)
    args = ["-ooc=%s.11.ooc"%fasta, "-dots=1000", "-noHead", "-extendThroughN", \
            "-minScore=%s"%identity, "-minIdentity=%s"%identity]
    cmd = "blat %s %s %s %s.psl"%(" ".join(args), fasta, fasta, fasta)
    if not verbose:
        cmd += " > /dev/null"
    else:
        sys.stderr.write(cmd+'\n')
    #generate overepresented 11mers if not exists
    if not os.path.isfile(fasta+".11.ooc"):
        os.system(cmd.replace("-ooc=", "-makeOoc="))
    #run BLAT
    os.system(cmd)
    #sort and take into account only larger vs smaller
    cmd2 = "awk '$10!=$14 && $11>=$15' %s.psl | sort -k11nr,11 -k12n,12 -k13nr,13 | gzip > %s.psl.gz"%(fasta, fasta)
    if verbose:
        sys.stderr.write(cmd2+'\n')
    os.system(cmd2)

def get_ranges(starts, sizes, offset=1):
    """Return str representation of alg ranges"""
    ranges = []
    for start, size in zip(starts.split(',')[:-1], sizes.split(',')[:-1]):
        start, size = int(start), int(size)
        start += offset
        end = start + size - offset
        coords = "%s-%s"%(start, end)
        ranges.append(coords)
    return " ".join(ranges)

def psl2hits(psl, identityTh, overlapTh, dblength=0, Lambda=0.318, K=0.13):
    """Return valid hits.
    PSL has to be sorted by q size: sort -k11nr,11 -k12n,12 -k13nr,13
    """
    hits = []
    added = set()
    for l in gzip.open(psl):
        if not l.strip() or not l.split()[0].isdigit():
            continue
        ##BLAT PSL without header
        (matches, mismatches, repm, Ns, Qgapc, Qgaps, Tgapc, Tgaps, strand, \
         q, qsize, qstart, qend, t, tsize, tstart, tend, blocks, bsizes, \
         qstarts, tstarts) = l.split('\t')
        #skip reverse matches - pairs tracking
        if q==t or (t,q) in added:
            continue
        added.add((q,t))
        #unpack batch
        matches, mismatches = int(matches), int(mismatches)
        Tgapc, Tgaps = int(Tgapc), int(Tgaps)
        qstart, qend = int(qstart), int(qend)
        qstart, qend, qsize = int(qstart), int(qend), int(qsize)
        tstart, tend, tsize = int(tstart), int(tend), int(tsize)
        #get score, identity & overlap
        score    = matches * 5 + mismatches * -3 + Tgapc * -4 + Tgaps * -1
        alglen   = int(tend) - int(tstart)
        identity = 1.0 * matches / alglen
        overlap  = 1.0 * alglen / tsize
        #filter by identity and overlap
        if overlap < overlapTh or identity < identityTh:
            continue
        #bitscore & evalue
        bitscore = (Lambda*score-math.log(K))/math.log(2)
        pvalue = evalue = 0
        if dblength:
            pvalue = 2**-bitscore
            evalue = len(qseqs[0]) * dblength * pvalue
        #get t and q ranges
        qranges  = get_ranges(qstarts, bsizes)
        tranges  = get_ranges(tstarts, bsizes)
        #store
        hits.append((q, qsize, qstart, qend, t, tsize, tstart, tend, \
                     identity, overlap, bitscore, evalue))
    return hits

def fasta2homozygous(out, fasta, identity, overlap, verbose):
    """Parse alignments and report homozygous contigs"""
    if verbose:
        sys.stderr.write("Indexing fasta...\n")
    #open fasta index
    faidx = SeqIO.index_db(fasta.name+".db3", fasta.name, "fasta")
    genomeSize = sum(len(faidx[c]) for c in faidx) 
    sys.stderr.write(" %s bp in %s contigs\n"%(genomeSize, len(faidx)))
    
    #run blat
    psl = fasta.name + ".psl.gz"
    if not os.path.isfile(psl):
        if verbose:
            sys.stderr.write("Running BLAT...\n")
        blat(fasta.name, identity, verbose)
    
    if verbose:
        sys.stderr.write("Parsing alignments...\n")
    #filter alignments
    matches = psl2hits(psl, identity, overlap)

    #remove redundant
    contig2skip = {}
    for c in faidx.keys():
        contig2skip[c] = 0
    for i, (q, qsize, qstart, qend, t, tsize, tstart, tend, identity, overlap, \
            bitscore, evalue) in enumerate(matches, 1):
        if q not in contig2skip:
            sys.stderr.write(' [ERROR] `%s` (%s) not in contigs!\n'%(q,str(matches[i-1])))
            continue
        #inform about matching already removed contig
        if verbose and contig2skip[q]:
            info = " [WARNING]: Match to already removed conting: %s %s\n"
            sys.stderr.write(info%(q, str(matches[i-1])))
        #store
        contig2skip[t] += 1

    #report homozygous fasta
    tsize = skipped = ssize = 0
    for i, c in enumerate(faidx, 1):
        if contig2skip[c]:
            skipped += 1
            ssize   += len(faidx[c])
            continue
        out.write(faidx[c].format('fasta'))
    #    
    sys.stderr.write(" Filtered out %s bp in %s heterozygous contigs!\n"%(ssize, skipped))
    sys.stderr.write("Reported %s bp in %s contigs.\n"%(genomeSize-ssize, len(faidx)-skipped))
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0a')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-f", "--fasta", required=True, type=file, 
                        help="FASTA file")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("--identity",    default=0.8, type=float, 
                        help="min. identity   [%(default)s]")
    parser.add_argument("--overlap",     default=0.75, type=float, 
                        help="min. overlap    [%(default)s]")
    #parser.add_argument("-p", "--ploidy",    default=2, type=int, 
    #                    help="ploidy          [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    #process fasta
    fasta2homozygous(o.output, o.fasta, o.identity, o.overlap, o.verbose)

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
