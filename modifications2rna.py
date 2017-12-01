#!/usr/bin/env python
desc="""Convert files with RNA modifications into canonical FASTA.
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Warsaw, 30/11/2017
"""

import sys
reload(sys)  
sys.setdefaultencoding('utf8')

from datetime import datetime
from collections import Counter

def table2modifications(table, log=sys.stderr):
    """Load modifications table"""
    if log: log.write("Loading modifications from %s...\n"%table)
    mod2base = {"U": "T", "_": "n", 'A': 'A', 'C': 'C', 'G': 'G'}
    mod2name = {}
    for l in open(table):
        #name	short_name	new_nomenclature	originating_base	rnamods_abbrev	html_abbrev	formula	monoisotopic_mass	average_mass
        ldata = l[:-1].split('\t')
        base, mod = ldata[3:5]
        if not l[:-1] or l.startswith('#') or not base or not mod or base==mod:
            continue
        mod2name[mod] = ldata[0]
        #if mod not in mod2base:    
        mod2base[mod] = base
    if log: log.write(" %s modifications loaded!\n"%len(mod2name))
    return mod2base, mod2name

def fasta_parser(fn):
    """Simple unicode aware FastA parser.
    Return header and seqs as a list
    """
    data = []
    for l in open(fn):
        l = l[:-1]
        if l.startswith(">"):
            if data:
                yield data
            data = [l[1:], []]
        else:
            data[-1] += [b.encode('utf-8') for b in l.upper().decode('utf-8')]
    if data:
        yield data
    
def get_id(desc, mod2base): 
    desc = desc.strip().replace(' | ', '.').replace(' ', '_')
    desc = "".join(b for b in desc if b.isalnum() or b in '._') # else '?' # or b not in mod2base else mod2base[b]
    return desc
    
def modifications2rna(table, fnames, species):
    """Convert RNA with modifications into canonical FastA"""
    mod2base, mod2name = table2modifications(table)

    print("Processing sequences...")
    II, MODS = 0, []
    i = 0
    for i, fn in enumerate(fnames, 1):
        ii = k = 0
        mods = []
        names = set()
        # define output files
        outfn1, outfn2 = fn+".rna.fa", fn+".dna.fa"
        print " %s --> %s & %s"%(fn, outfn1, outfn2)
        with open(outfn1, "w") as out1, open(outfn2, "w") as out2:
            for name, seq in fasta_parser(fn):
                # skip if not given species
                if species and species not in name:
                    continue
                ii += 1
                # get name
                name = get_id(name, mod2base)
                while name in names: name+="_"
                names.add(name)
                # remove gaps
                seq = filter(lambda x: x!="-", seq)
                seq1 = "".join(seq)
                # count modifications
                mods += filter(lambda b: b in mod2name, seq)
                # recode to DNA # replace is needed as some modifications are replaced to U
                seq2 = "".join(mod2base[b] for b in seq).replace('U', 'T') 
                if len(seq1.decode('utf-8'))!=len(seq2.decode('utf-8')):
                    k += 1
                    print ii, k, name; print seq1; print seq2#; print [mod2base[b] for b in seq]
                out1.write(">%s\n%s\n"%(name, seq1))
                out2.write(">%s\n%s\n"%(name, seq2))
        # report mod2count
        report_stats(fn, ii, mods, mod2name, species)
        II += ii
        MODS += mods
    # report summary only if more than 1 file
    if i>1:
        report_stats("all", II, MODS, mod2name, species)
        
def report_stats(fn, ii, mods, mod2name, species="", log=sys.stderr):
    """Report stats about modifications"""
    speciesinfo = ""
    if species:
        speciesinfo = " for %s"%species
    mod2count = Counter(mods)
    log.write("# %s has %s modifications of %s types in %s seqs%s:\n" % (fn, len(mods), len(mod2count), ii, speciesinfo))
    for mod, c in mod2count.most_common(50):
        log.write("%s\t%s\t%s\n"%(mod, c, mod2name[mod]))

def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.1') 
    parser.add_argument("-i", "--fnames", nargs="+", help="files to process" )
    parser.add_argument("-t", "--table", default="modifications.txt", help="modification table [%(default)s]" )
    parser.add_argument("-s", "--species", default="", help="limit entries to those having given str [%(default)s]" )
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write( "Options: %s\n" % str(o) )

    modifications2rna(o.table, o.fnames, o.species)

if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
    