#!/usr/bin/env python
desc="""Insert SNPs from vcf into fasta.
If coordinates are provided in header name (chrI:200-300) the script will use them!
NOTE: coordinates are 0-based, half-open (BED,python) by default. 
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Girona/Cracow, 6/11/2012
"""

import argparse, os, random, sys
from datetime import datetime
from Bio      import SeqIO

def load_vcf(vcf, noindels, nohetero, verbose):
    """Return dictionary of variants
    with 0-based, half-open coordinates (BED,python).
    VCF is 1-based (GFF like) ?
    """
    chr2variants = {}
    for l in vcf:
        l = l.strip()
        if not l or l.startswith("#") or "PASS" not in l:
            continue
        contig, pos, id, ref, alt, qual, flt = l.split("\t")[:7]
        #0-based coordinates
        pos = int(pos)-1
        if noindels and len(ref)!=len(alt):
            continue
        #skip heterozygous SNPs if defined so
        if '0/1' in l:
            #or select random
            if nohetero or random.randint(0, 1):
                continue
        if not set(alt).issubset(set('ACGT')):
            info = " Warning: weird alternative base: %s:%s %s %s\n"
            sys.stderr.write(info%(contig,pos,ref,alt))
            continue
        #add chr to dict
        if not contig in chr2variants:
            chr2variants[contig]=[]
        #append
        chr2variants[contig].append((pos, ref, alt))
    return chr2variants

def variants2sequence(gene, seq, vof, qsorg, verbose):
    """Introduce variants into sequence.
    Variants positions has to be 0-based (BED,python).
    TO DO: here don't use pos but qsorg - pos
    """
    #introduce variants from the end, so indels will not mess up coordinates
    snps = indels = 0
    #bigindels = []
    for opos, ref, alt in sorted(vof, reverse=True):
        #update pos
        pos = opos - qsorg
        #warn about indel
        if len(ref)!=len(alt):
            indels += 1
        else:
            snps   += 1
        #check ref
        if str(seq[pos:pos+len(ref)]) != ref:
            info = " Warning: variants2sequence: Orginal seq (%s) differs from variant ref (%s) in %s!\n"
            sys.stderr.write(info % (seq[pos:pos+len(ref)], ref, gene))
        #replace orginal seq
        #print seq
        seq = seq[:pos] + alt + seq[pos+len(ref):]
        if verbose:
            sys.stderr.write( "%s: %s > %s\n" % (pos, ref, alt))
    return (str(seq), snps, indels)
    
def vcf2fasta(out, vcf, genome, noindels, nohetero, verbose):
    """
    """
    #load vcf file
    chr2variants = load_vcf(vcf, noindels, nohetero, verbose)

    #process all entries
    for r in SeqIO.parse(genome, "fasta"):
        q = r.id
        vof = []
        #carefully check if chr1:100-200 or just chr1
        if len(q.split(":"))==2 and len(q.split(":")[1].split("-"))==2 \
           and q.split(":")[1].split("-")[0].isdigit() \
           and q.split(":")[1].split("-")[1].isdigit():
            #get original seq start end end
            qcontig     = q.split(":")[0]
            qsorg,qeorg = q.split(":")[1].split("-")
            qsorg,qeorg = int(qsorg),int(qeorg)
            #get overlapping variants
            if qcontig in chr2variants:
                vof = filter(lambda x: qsorg <= x[0] < qeorg,chr2variants[qcontig])
        else:
            #define start pos
            qsorg = 0
            if q in chr2variants:
                vof = chr2variants[q]
        #introduce variants into seq    
        seq, snps, indels = variants2sequence(r.id, r.seq, vof, qsorg, verbose)
        #report big difference in seq length after snp/indel correction
        if abs(len(seq)-len(r.seq)) > 5:
            info = " Warning: Gene %s: %s difference between orginal (%s) and corrected (%s) sequence length!\n"
            sys.stderr.write(info % (r.id, abs(len(seq)-len(r.seq)), len(seq), \
                                     len(r.seq)))
        #report entry with variants
        line = ">%s snps:%s indels:%s\n%s\n"
        out.write(line % (r.id, snps, indels, \
                          "\n".join(seq[i:i+60] for i in range(0, len(seq), 60))))
        
def main():

    usage  = "usage: %(prog)s [options]"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
    
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="vcf",     default=sys.stdin, type=file, 
                        help="vcf file with SNPs  [stdin]")
    parser.add_argument("--noindels",        default=False, action="store_true", 
                        help="skips INDELs")    
    parser.add_argument("--nohetero",        default=False, action="store_true", 
                        help="skip heterozygous SNPs")    
    parser.add_argument("-g", dest="genome",  type=file, 
                        help="genome fasta file   [mandatory]")
    parser.add_argument("-o", dest="out",     default=sys.stdout, type=argparse.FileType('w'), 
                        help="define output name  [stdout]")
    #parser.add_argument("-1", dest="onebased",default=False, action="store_true", 
    #                    help="1-based coordinates [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    #go
    vcf2fasta(o.out, o.vcf, o.genome, o.noindels, o.nohetero, o.verbose)


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
