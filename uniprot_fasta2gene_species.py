#!/usr/bin/env python3
# Extract protid, gene_name and species_name from Uniprot FastA
# USAGE: cat uniprot.fasta | uniprot_fasta2gene_species.py

import sys

for l in sys.stdin:
    if l.startswith(">"):
        #>tr|A0A2K5WXK9|A0A2K5WXK9_MACFA Uncharacterized protein OS=Macaca fascicularis OX=9541 GN=RPL3L PE=3 SV=1 > A0A2K5WXK9 RPL3L Macaca_fascicularis
        protid = l.split('|')[1]
        k2v = {} #: 
        k = ""
        for kv in l.split("="):
            if not k:
                k = kv.split()[-1]
            else:
                v = " ".join(kv.split()[:-1]).strip()
                k2v[k] = v
                k = kv.split()[-1]
        species_name = k2v["OS"]
        gene_name = k2v["GN"]
        print("\t".join([protid, species_name, gene_name]))
        
        
