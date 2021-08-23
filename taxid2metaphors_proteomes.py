#!/usr/bin/env python2
desc="""Fetch proteomes of taxa from given group form metaPhOrs. 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Mizerow, 17/07/2014
"""

import argparse, os, sys
from datetime import datetime
sys.path.insert(0, '/users/tg/lpryszcz/cluster/rapsi/src')
from taxonomy import Taxonomy
sys.path.insert(0, '/users/tg/lpryszcz/cluster/metaphors/src/client')
import dbClient

def taxid2proteomes(m, taxa, group_taxid, verbose):
    """Fetch proteomes for given taxid"""
    #get group info
    parent_id, name, rank = taxa.get_taxa_info(group_taxid)
    outdir = "%s_%s"%(group_taxid, name) #name
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    #
    taxids = filter(lambda taxid: group_taxid in taxa.taxid2lineage(taxid), m.species)
    #taxids = filter(lambda taxid: name in m.species[taxid][-1], m.species)
    if verbose:
        info = "Fetching %s taxa proteomes for %s [%s %s]...\n"
        sys.stderr.write(info%(len(taxids), group_taxid, name, rank))
        #sys.stderr.write("Fetching %s taxa for %s...\n"%(len(taxids), name))
    #save fasta
    for i, taxid in enumerate(taxids, 1):
        spname = m.species[taxid][1]
        sys.stderr.write(" %s  %s                       \r"%(i, spname))
        outfn = os.path.join(outdir, "%s.%s.faa"%(taxid, spname.replace(' ','_')))
        with open(outfn, 'w') as out:
            out.write(m.get_proteome(taxid))
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")
    parser.add_argument('-d', '--db', default="metaphors_201405",
                        help="database name     [%(default)s]")
    parser.add_argument('-t', '--taxids', nargs="+", type=int,
                        help="group taxid(s)    [%(default)s]")
    parser.add_argument("--taxadb",        default="/users/tg/lpryszcz/cluster/rapsi/taxonomy.db3",
                        help="taxonomy path  [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    #init taxonomy
    taxa = Taxonomy(o.taxadb)

    #init metaphors connection
    m  = dbClient.metaphors(o.db)
    if o.verbose:
        sys.stderr.write("%s species in %s database\n"%(len(m.species), o.db))
        
    #process taxa groups
    for taxid in o.taxids:
        #fetch proteins from given taxa
        taxid2proteomes(m, taxa, taxid, o.verbose)
        
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
