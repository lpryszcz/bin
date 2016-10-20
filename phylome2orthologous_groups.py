#!/usr/bin/env python
desc="""Return orthologous groups for given species from a phylome.
It adds orthologs (and in-paralogs) to the group until outgroup species found.

TBD: Add one2one testing!
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 8/07/2013
"""

import argparse, os, sys
from datetime import datetime
from ete2     import PhyloTree, PhylomeDB3Connector
from phylome  import _getConnection, _get_spcode, get_evolEvents
from rooted_phylomes import ROOTED_PHYLOMES

def root_tree(t, seedid, phyid):
    """Return rooted tree"""
    #get seednode
    seedNode = t.get_leaves_by_name(seedid)[0]
    #ROOT TREE with error avoinding loops - should be sorted out in the future!    
    try:
        if phyid in ROOTED_PHYLOMES:
            t.set_outgroup(seedNode.get_farthest_oldest_node(ROOTED_PHYLOMES[phyid]))
        else:
            t.set_outgroup(t.get_midpoint_outgroup())
    except:
        sys.stderr.write( "  Warning: get_evolEvents: Cannot root tree: %s %s\n" % (seedid,t.write()) )
        return None
    return t
    
def tree2orthogroups(t, seedid, species, one2one, fraction, verbose):
    """Return orthogroup for given species from tree
    Add one2one testing!"""
    #find biggest subpartition containing seed and only species from species
    ogroup, ospecies = [], []
    for i, subtree in enumerate(t.traverse('postorder'), 1):\
        #skip subpartitions without seed
        if seedid not in subtree:
            continue
        #break if any species in partition from outside species
        if subtree.get_species().difference(species):
            break
        #store orthogroup
        ogroup   = subtree.get_leaf_names()
        ospecies = subtree.get_species()
    #check if fraction passed
    if len(ospecies)*1.0/len(species) < fraction:
        return None
    return sorted(ogroup)

def phylome2orthogroups(out, phyid, species, one2one, fraction, step, verbose):
    """Report orthogroups for given phylome and set of species"""
    #get phylomeDB connection
    p = _getConnection()
    #get seeds
    seedids = p.get_phylome_seed_ids(phyid)[0]
    if verbose:
        sys.stderr.write("Processing trees...\n")
    processed = set()
    trees = rooted = ogroups = 0
    for i, seedid in enumerate(seedids, 1):
        if seedid in processed:
            continue
        #get tree
        trees_dict = p.get_best_tree(seedid, phyid)
        if not trees_dict:
            continue
        t = trees_dict['tree']
        if not t:
            continue
        trees += 1
        #set species naming function
        t.set_species_naming_function(_get_spcode)         
        #root
        #t = root_tree(t, seedid, phyid)
        if not t:
            continue
        rooted += 1
        #get orthogroup and report
        orthogroup = tree2orthogroups(t, seedid, species, one2one, fraction, verbose)
        if not orthogroup:
            continue
        ogroups += 1
        out.write("\t".join(orthogroup)+"\n")
        #add to processed
        for o in orthogroup:
            processed.add(o)
        #print info
        if verbose and i%step == 1:
            sys.stderr.write(" %s / %s %s %s %s\r" % (i, len(seedids), trees, rooted, ogroups))

def main():
    usage   = "%(prog)s -v" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",   default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-p", "--phyid",       default=0, type=int,
                        help="define phylome id                 [%(default)s]")
    parser.add_argument("-s", "--species",      nargs="+", 
                        help="list of species")
    parser.add_argument("-o", "--out",          default=sys.stdout, type=argparse.FileType("w"), 
                        help="define output stream              [stdout]")
    parser.add_argument("-f", "--fraction",     default=0.0, type=float,
                        help="fraction of species in orthogroup [%(default)s]")
    parser.add_argument("--one2one",            default=False, action="store_true",
                        help="allowed only one2one orthologs (not implemented)")
    parser.add_argument("--step",               default=100, type=int,
                        help="print info every t trees          [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stdout.write( "[gi2gene] Options: %s\n" % str(o) )

    phylome2orthogroups(o.out, o.phyid, o.species, o.one2one, o.fraction, o.step, o.verbose)
    
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\n Ctrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stdout.write( "#Time elapsed: %s\n" % dt )
