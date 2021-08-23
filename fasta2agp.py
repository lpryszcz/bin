#!/usr/bin/env python2
desc="""Generate AGP from scaffolds FASTA.
"""
epilog="""AUTHOR:
l.p.pryszcz+git@gmail.com

Mizerow, 27/02/2014
"""

import argparse, os, re, sys
from datetime import datetime
from Bio import SeqIO

gappat = re.compile('N+')
def fasta2agp(handle, outbase, minN, evidence, verbose):
    """Generate AGP 2.0
    https://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/AGP_Specification.shtml
    """
    #define gap pattern
    pat = re.compile('N{%s,}'%minN)
    #open out files
    contigs = open(outbase+".contigs.fa", "w")
    agp     = open(outbase+".agp", "w")
    #parse fasta
    line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
    for i, r in enumerate(SeqIO.parse(handle, 'fasta'), 1):
        if verbose:
            sys.stderr.write(" %s %s %s bp      \r"%(i, r.id, len(r)))
        pend = 0
        j = 1
        for m in pat.finditer(str(r.seq).upper()):
            #get gap span
            gstart, gend = m.span()
            #write contig
            c = r[pend:gstart]
            c.id = "%s.contig%s"%(r.id, j)
            c.description = c.id
            contigs.write(c.format("fasta"))
            #write agp
            #agp.write(line%(r.id, pend+1, gstart, 2*j-1, "W", c.id, pend+1, gstart, "+"))
            agp.write(line%(r.id, pend+1, gstart, 2*j-1, "W", c.id, 1, gstart-pend, "+"))
            agp.write(line%(r.id, gstart+1, gend, 2*j, "N", gend-gstart, "scaffold", "yes", evidence))
            #store previous end
            pend = gend
            j += 1
        #store last bit
        gstart = len(r)
        #write contig
        c = r[pend:gstart]
        c.id = "%s.contig%s"%(r.id, j)
        c.description = c.id 
        contigs.write(c.format("fasta"))
        #write agp
        #agp.write(line%(r.id, pend+1, gstart, 2*j-1, "W", c.id, pend+1, gstart, "+"))
        agp.write(line%(r.id, pend+1, gstart, 2*j-1, "W", c.id, 1, gstart-pend, "+"))
    #close output streams
    agp.close()
    contigs.close()

def main():
    usage   = "%(prog)s -o output_basename"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
    
    parser.add_argument("-v", "--verbose",  default=False, action="store_true")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument("-i", "--input", default=sys.stdin, type=file, 
                        help="input stream  [stdin]")
    parser.add_argument("-o", "--output_base", required=True, #default=sys.stdout, type=file,
                        help="output base name")
    parser.add_argument("-n", "--minN", default=10, type=int, 
                        help="min N to call gap [%(default)s]")
    parser.add_argument("-e", "--evidence", default="paired-ends", 
                        help="evidence type [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n" % str(o))

    fasta2agp(o.input, o.output_base, o.minN, o.evidence, o.verbose)

if __name__=='__main__': 
    t0=datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt=datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
