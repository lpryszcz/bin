#!/usr/bin/env python
desc="""Detect premature termination codons and frame-shifts from exonerate alignments.

Run exonerate with this parameteres
    exonerate --model protein2genome -n 1 --softmaskquery no --softmasktarget yes --minintron 20 --maxintron 20000 --geneseed 250 --query Ir.pep.fa --target ../../ref/$s.fa > Ir.pep.fa.$s.exout1

"""
epilog="""Author: l.p.pryszcz+github@gmail.com
Carmona, 29/03/2017
"""

import os, sys
from datetime import datetime
from FastaIndex import FastaIndex

def _get_alignement(handle):
    """Return alignement"""

def parse_exonerate(handle):
    """parse exonerate output"""
    alg = ['', '', '', '']
    ii = 0
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
            ii = 0
            ldata = l.split(': ')[-1].split()
            q, qs, qe, qstrand, t, ts, te, tstrand = ldata[:8]#; print ldata
            ts, te, qs, qe = map(int, (ts, te, qs, qe))
            yield q, qs, qe, qstrand, t, ts, te, tstrand, alg
        # read alignment
        else:
            if ii%4 in (0, 3):
                l = l.split(' : ')[1]
                pl = l    
            else:
                l = l[-len(pl):]
            alg[ii%4] += l
            ii += 1    
            

def exonerate2ptc(handle, out, pepfn, minoverlap, verbose):
    """Return premature stop codons and frame-shifts from exonerate output"""
    faidx = FastaIndex(pepfn)
    pep2len = {pep: data[0] for pep, data in faidx.id2stats.iteritems()}
    for q, qs, qe, qstrand, t, ts, te, tstrand, alg in parse_exonerate(handle):
        alglen =  map(len, alg)
        if len(set(alglen))>1:
            print "error", q, t, alglen
        if (qe-qs)/pep2len[q] < minoverlap:
            print t, q, alglen#; print '\n'.join(alg); return
            continue
        

def main():
    import argparse
    usage = "%(prog)s -i " #usage=usage, 
    parser = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                          formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--version', action='version', version='0.11c')	 
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")	
    parser.add_argument("-i", "--input", default=sys.stdin, type=file, help="input stream [stdin]")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream [stdout]")
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
    