#!/usr/bin/env python
desc="""Detect premature termination codons and frame-shifts from exonerate alignments.

Run exonerate with this parameteres
 exonerate --model protein2genome -n 1 --softmaskquery no --softmasktarget yes --minintron 20 --maxintron 20000 --query Ir.pep.fa --target ../../ref/$s.fa > Ir.pep.fa.$s.exout1

"""
epilog="""Author: l.p.pryszcz+github@gmail.com
Carmona/Malaga/Brussels/Warsaw, 29/03/2017
"""

import os, sys, re
from datetime import datetime
from FastaIndex import FastaIndex
            
ptcpat = re.compile('\*{3}')
frameshiftpat = re.compile('#+')
intronpat = re.compile('\s+>>>> Target Intron \d+ >>>>\s+')
    
def parse_exonerate(handle):
    """parse exonerate output"""
    ii = 0
    alg = ['', '', '', '']
    symbols = ("-", " ", " ", ".")
    while True:
        l = handle.readline()
        if not l: break
        l = l[:-1]#; print l
        # skip header and spacer lines
        if l.startswith(('Command line:', 'Hostname:', '--')) or not l:
            continue
        # alg start - skip header
        if l.startswith('C4 Alignment:'):
            while l:
                l = handle.readline()[:-1]
        # read vulgar and report
        elif l.startswith('vulgar: '):
            ldata = l.split(': ')[-1].split()
            q, qs, qe, qstrand, t, ts, te, tstrand, score = ldata[:9]#; print ldata
            vulgar = ldata[9:]
            ts, te, qs, qe = map(int, (ts, te, qs, qe))
            # remove {}
            for i in range(len(alg)):
                alg[i] = alg[i].replace('{', '').replace('}', '')
            # get ipos and isizes and 
            ipos = reversed([m.span() for m in intronpat.finditer(alg[0])])
            isizes = reversed([int(qi)+4 for v, qi in zip(vulgar[::3], vulgar[2::3]) if v=="I"])
            for isize, (istart, iend) in zip(isizes, ipos):
                for i in range(len(alg)):
                    alg[i] = alg[i][:istart] + symbols[i]*isize + alg[i][iend:]
            # recode <> to -
            for i in (0, 3):
                alg[i] = alg[i].replace('<', '-').replace('>', '-')
            # yield
            yield q, qs, qe, t, ts, te, tstrand, score, vulgar, alg
            # reset
            ii = 0
            alg = ['', '', '', '']
        # read alignment
        else:
            if ii%4 in (0, 3):
                l = l.split(' : ')[1]
                pl = l    
            else:
                l = l[-len(pl):]
            alg[ii%4] += l
            ii += 1

def get_position(seq, sstart, send, strand, s, e, posdivide=1):
    """Return feature position"""
    if strand in "+":
        pos = sstart + s - seq[:s].count("-")
    else:
        pos = send*posdivide + sstart-send - e - seq[:s].count("-")
    return pos/posdivide
    
def print_alg(q, t, alg, linelen=120):
    """Print formatted alignment"""
    print q, t
    for s in range(0, len(alg[0]), linelen):
        for i in range(len(alg)):
            print alg[i][s:s+linelen]
        print
    
def exonerate2ptc(handle, out, pepfn, minoverlap, verbose):
    """Return premature stop codons and frame-shifts from exonerate output"""
    # get peptide lengths
    faidx = FastaIndex(pepfn)
    pep2len = {pep: data[0] for pep, data in faidx.id2stats.iteritems()}
    # process exonerate matches
    i = k = 0
    ptc, frameshift = [], []
    lowoverlap = " %s - %s:%s-%s%s skipped due to low overlap (%.3f)\n"
    out.write("# chrom\tstart\tend\tgene aaPOSalt\tscore\tstrand\n")
    for i, (q, qs, qe, t, ts, te, strand, score, vulgar, alg) in enumerate(parse_exonerate(handle), 1):
        gene_name = q.split('|')[-1]
        alglen =  map(len, alg)
        if len(set(alglen))>1:
            sys.stderr.write("[ERROR] Uneven alignments lengths: %s - %s %s\n" %(q, t, str(alglen)))
            continue
        # check overlap
        overlap = 1.*(qe-qs)/pep2len[q]
        if overlap < minoverlap:
            if verbose: sys.stderr.write(lowoverlap%(q, t, ts, te, tstrand, overlap))
            continue
        # unload alg
        qaminos, msymbols, taminos, tseq = alg
        # report ptc
        for ii, m in enumerate(ptcpat.finditer(taminos[:-3]), 1):
            ptc.append(q)
            s, e = m.span()
            # get amino positon - always +
            qpos = get_position(alg[0], qs, qe, "+", s, e, 3)
            name = "%s %s%s*" % (gene_name, alg[0][s:e], qpos+1)
            # get chrom position
            tpos = get_position(alg[-1], ts, te, strand, s, e, 1)
            bed = "%s\t%s\t%s\t%s\t%.3f\t%s\n"%(t, tpos, tpos+e-s, name, 1-1.*qpos/pep2len[q], strand)
            out.write(bed); 
        # report frame shifts
        for ii, m in enumerate(frameshiftpat.finditer(taminos), 1):
            frameshift.append(q)
            s, e = m.span()
            qpos = get_position(alg[0], qs, qe, "+", s, e, 3)
            name = "%s %sdel%s" % (gene_name, qpos+1, e-s) #alg[0][s:e], 
            # get chrom position
            tpos = get_position(alg[-1], ts, te, strand, s, e, 1)
            bed = "%s\t%s\t%s\t%s\t%.3f\t%s\n"%(t, tpos, tpos+e-s, name, 1-1.*qpos/pep2len[q], strand)
            out.write(bed)
        #if gene_name=="Ir60d": print_alg(q, t, alg)
        k += 1

    # get unique
    ptcset, frameshiftset = set(ptc), set(frameshift)
    sys.stderr.write("Processed %s alignments for %s proteins, of these %s had overlap >= %.3f\n"%(i, len(pep2len), k, minoverlap))
    sys.stderr.write(" %s PTCs in %s proteins: %s\n"%(len(ptc), len(ptcset), " ".join(ptcset)))
    sys.stderr.write(" %s frameshifts in %s proteins: %s\n"%(len(frameshift), len(frameshiftset), " ".join(frameshiftset)))

def main():
    import argparse
    usage = "%(prog)s -i " #usage=usage, 
    parser = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                          formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--version', action='version', version='0.11c')	 
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")	
    parser.add_argument("-i", "--input", default=sys.stdin, type=file, help="input stream [stdin]")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType('w'), help="output stream [stdout]")
    parser.add_argument("-p", "--pep", required=1, help="peptides FastA file")
    parser.add_argument("--overlap", default=0.9, help="min. fraction of query aligned [%(default)s]")

    '''# print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)'''
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    exonerate2ptc(o.input, o.output, o.pep, o.overlap, o.verbose)
	
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!		\n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    #[Errno 95] Operation not supported
    except OSError:
        sys.stderr.write("OS error({0}): {1}\n".format(e.errno, e.strerror))        
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
    