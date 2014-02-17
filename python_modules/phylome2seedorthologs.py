#!/usr/bin/env python
desc="""Process trees from phylome tar.gz dump
and generates orthologs dump of seed species.

KNOWN ISSUES:
- if one tree is wrongly rooted, co-orthologs are shared by all trees
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 11/01/2013
"""

import argparse, os, sys
import gzip, tarfile
import numpy, resource
from datetime import datetime
from ete2     import PhyloTree
from rooted_phylomes  import ROOTED_PHYLOMES

def phylomedump_tree_iterator( tarfn,verbose=0 ):
    """PhylomeDB all_trees.tar.gz dump treeobj generator."""
    #open tarfile
    if tarfn.endswith(".gz"):
        tar = tarfile.open(tarfn,"r:gz")
    else:
        tar = tarfile.open(tarfn,"r")

    i = k = 0
    #process entries
    for m in tar:
        #if i>100: break
        if not m.isfile():
            continue
        #load tree
        if   m.name.endswith(".nw"):
            i += 1
            #get nw
            nw = tar.extractfile(m).readline()
            t  = PhyloTree(nw)
            ##add seedid and method info
            #Phy000CWA9_YEAST.JTT.nw --> Phy000CWA9_YEAST JTT
            seedid, method = os.path.basename(m.name).split(".")[:2]
            t.seedid = seedid
            t.method = method
        #or add lk, seedid, method and lk to treeobj
        elif m.name.endswith(".lk"):
            seedid, method, lk = tar.extractfile(m).readline().split('\t')[:3]
            t.lk = float(lk)
            if not t.lk:
                sys.stderr.write( " Err: Zero likelihood (%s) for: %s\n" % (t.lk, ", ".join((t.seedid, t.method))))
                continue
            if seedid!=t.seedid or t.method != method:
                sys.stderr.write( " Err: Seedid and/or method doesn't match: %s\n" % ", ".join((seedid, t.seedid, method, t.method)))
                continue
            k += 1
            if verbose and not i%100:
                sys.stderr.write( "  %6i\r" % i )

            yield t
    if verbose:
        sys.stderr.write( " %s out of %s trees succesfully parsed [memory: %s KB]\n" % (k, i, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

def _get_spcode( leaf ):
    """Return species code"""
    code = leaf.split('_')[-1]
    return code

def add_homologs(homologs, seedspecies, lk, h1list, h2list, ee, verbose):
    """
    """
    for h1 in h1list:
        #skip non-seedsp proteins
        if _get_spcode(h1) != seedspecies:
            continue
        for h2 in h2list:
            #skip within-species paralogs
            if _get_spcode(h2) == seedspecies:
                continue
            #prot1,prot2 = h1,h2
            pair = "%s-%s" % (h1, h2)
            #print prot1,prot2,lk
            if pair not in homologs:
                homologs[pair] = []
            homologs[pair].append((lk, ee))
    return homologs

def tree2seedhomologs(t, homologs, seedspecies, lk, verbose):
    """Add homologs from given tree to dict"""
    evolEvents = t.get_descendant_evol_events()
    #GET EVOLEVENTS
    for e in evolEvents:
        #define species overlap
        ee = 1
        if not e.sos:
            ee = 0
        #add seedsp homologs
        homologs = add_homologs(homologs, seedspecies, lk, e.in_seqs, e.out_seqs, ee, verbose)
        homologs = add_homologs(homologs, seedspecies, lk, e.out_seqs, e.in_seqs, ee, verbose)
    return homologs

def process_trees(tarfn, phyid, verbose):
    """Add homologies from given tree to dict."""
    homologs = {}
    for t in phylomedump_tree_iterator(tarfn, verbose):
        seedid = t.seedid
        species = _get_spcode(seedid)
        method = t.method
        lk     = t.lk
        treeid = "%s_%s_%s" % (phyid, seedid, method)
        #get seedNode
        try:
            seedNode = t.get_leaves_by_name(seedid)[0]
        except:
            sys.stderr.write(" Err: Cannot get seedid leaf in: %s\n" % treeid)
            continue

        #ROOT TREE with error avoinding loops - should be sorted out in the future!
        t.set_species_naming_function(_get_spcode)
        if phyid in ROOTED_PHYLOMES:
            try:
                outgroup = seedNode.get_farthest_oldest_leaf(ROOTED_PHYLOMES[phyid])
                t.set_outgroup(outgroup)
            except:
                sys.stderr.write(" Err: Cannot root tree: %s using sp2age\n" % treeid)
                continue
        else:
            try:
                t.set_outgroup(t.get_midpoint_outgroup())
            except:
                sys.stderr.write(" Err: Cannot root tree: %s using mid-point\n" % treeid)
                continue

        #get homologs of seed sequence
        homologs = tree2seedhomologs(t, homologs, species, lk, verbose)
        #if seedid == "Phy0064KBX_420245": print treeid; t.show()
    if verbose:
        sys.stderr.write(" %s homolog pairs stored\n" % len(homologs))
    return homologs

def homologs2orthologs(homologs, CSth, lkTh, verbose):
    """Return orthologs"""
    orthologs   = {}
    coorthologs = {}
    for pair in homologs.keys():
        #unpack pair
        pairList = homologs.pop(pair)
        p1,p2 = pair.split("-")

        #get bestLK for given pair
        bestLk = max([lk for lk, ee in pairList])
        #skip trees having lk < 3x bestLk
        ees = [ee for lk, ee in pairList if lk >= lkTh * bestLk]
        #check CS
        CS = 1 - numpy.mean(ees)
        if CS < CSth:
            continue

        #add ortholog entry
        sp2 = _get_spcode(p2)
        p1sp2 = '%s-%s' % (p1,sp2)
        #add seedid and sp to orthologs
        if p1sp2 not in orthologs:
            orthologs[p1sp2] = []

        #store orthologs data
        orthologs[p1sp2].append((p2, CS, len(ees)))

        #add co-orthologs entry
        if p2 not in coorthologs:
            coorthologs[p2] = []
        #store co-orthologs
        coorthologs[p2].append(p1)

    return orthologs, coorthologs

def orthologs2outfile(orthologs, coorthologs, outfn, verbose):
    """Store orthologs to out file."""
    header='#seqid\torthologid\ttype\tCS\ttrees\tco-orthologs\n'
    out   = gzip.open(outfn, 'w')
    out.write(header)
    for p1sp2 in orthologs:
        p1, sp2 = p1sp2.split('-')
        #get relation2 type
        rel2 = 'one'
        if len(orthologs[p1sp2])>1:
            rel2 = 'many'
        #parse orthologs from sp2
        for p2, CS, tree in orthologs[p1sp2]:
            co_orthologs_text = " ".join([co for co in coorthologs[p2] if co != p1]) #co_orthologs_text[:-1]
            #and get relation1 type
            rel1 = 'one'
            if co_orthologs_text:
                rel1 = 'many'
            #save line
            line='%s\t%s\t%s-to-%s\t%1.3f\t%i\t%s\n' % (p1, p2, rel1, rel2, CS, tree, co_orthologs_text)
            out.write(line)
    out.close()

def phylomedump2seedhomologs(tarfn, outfn, phyid, CSth, lkTh, verbose):
    """ """
    if verbose:
        sys.stderr.write("[%s] Parsing trees from: %s\n" % (datetime.ctime(datetime.now()), tarfn))
    homologs = process_trees( tarfn,phyid,verbose )

    #process all homologies and get co-homologies
    if verbose:
        sys.stderr.write("[%s] Computing meta predictions...\n" % datetime.ctime(datetime.now()) )
    orthologs,coorthologs = homologs2orthologs(homologs, CSth, lkTh, verbose)

    #store to file
    if verbose:
        sys.stderr.write( "[%s] Saving to outfile: %s [memory: %s KB]\n" % (datetime.ctime(datetime.now()), outfn, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    orthologs2outfile( orthologs,coorthologs,outfn,verbose )

def main():
    usage   = "%(prog)s [options] -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)

    parser.add_argument("-v", dest="verbose",   default=False, action="store_true", help="verbose")
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="infn",      required=True,
                        help="input file name   ")
    parser.add_argument("-p", dest="phyid",     required=True, type=int,
                        help="phylome id        [%(default)s]")
    parser.add_argument("-o", dest="outfn",     default="orthologs.txt.gz",
                        help="output file name  [%(default)s]")
    parser.add_argument("-c", dest="CS",        default=0.5, type=float,
                        help="CS cutoff         [%(default)s]")
    parser.add_argument("-l", dest="lk",        default=3.0, type=float,
                        help="likelihood cutoff [%(default)s]" )
    #Salva's patch
    parser.add_argument("-r", dest="rooting_dict", default="", type=str,
                        help="input file containing rooting dictionary" )
    #end
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n" % str(o))

    #Salva's patch
    global ROOTED_PHYLOMES

    if os.path.isfile(o.rooting_dict):
        rooting_dict = eval(open(o.rooting_dict, "rU").read())
        #for k in rooting_dict:
        #    ROOTED_PHYLOMES[k] = rooting_dict[k]
        ROOTED_PHYLOMES.update(rooting_dict)
    #end
    
    phylomedump2seedhomologs(o.infn, o.outfn, o.phyid, o.CS, o.lk, o.verbose)

if __name__=='__main__':
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
