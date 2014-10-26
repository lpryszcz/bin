#!/usr/bin/env python 
desc="""Get phylogenetic profile of orthologs for given list of proteins.

problems to resolve:
-if getting tree by seedfriend, you have to check if seed is ortholog of query sequence!

Code based on neptune:~/assembly_trial/virulence/profile.py
"""
epilog="""Author: l.p.pryszcz@gmail.com
Barcelona/Mizerow, 14/02/2012
"""

import ete2, gzip, os, sys
#sys.path.append('/home/lpryszcz/workspace/phylomeDB/src/') #sys.path.append('/home/lpryszcz/workspace/phylomeDB/src/')
#import phylomeDB3
from rooted_phylomes import ROOTED_PHYLOMES
from datetime import datetime
import numpy as np
from Bio import SeqIO
#from MyGraph import MyGraph

from phylome import _get_spcode, _getConnection

def get_protid2phyid(protid2phylmeFn, protids, spCode, p, phyid):
    """Load protd2phyid from file or from phylomeDB"""
    protid2phyid={}
    #from file if specified
    if protid2phylmeFn:
        for line in open(protid2phylmeFn):
            protid,phyprot=line[:-1].split()[:2]
            protid2phyid[protid]=phyprot
        return protid2phyid

    #get version of proteome used in phylome
    version = int(filter(lambda x: x[1]==spCode, p.get_proteomes_in_phylome(phyid))[0][2])
        
    #otherwise try to read from phylomedb
    k = 0
    info = "[WARNING] Cannot get phylome protid for: %s ( p.get_id_by_external(%s) -> %s )\n"
    for i, protid in enumerate(protids, 1):
        #get id
        '''#geneid=protid.split('.')[0]
        idData=p.get_id_by_external(protid) #{'CANPA': {1: ['Phy0022W7C_CANPA']}, 'COCPO': {1: ['Phy002ZINX_COCPO']}}
        phyprot=None
        if spCode in idData: 
            phyprot=idData[spCode].itervalues().next()[0]
        '''
        idData = p.search_id(protid)
        if version in idData:
            phyprot = idData[version]
        if not phyprot:
            sys.stderr.write(info%(protid, protid, idData))
            continue
        k += 1
        protid2phyid[protid] = phyprot
    sys.stderr.write("Found %s mappings for %s proteins.\n"%(k, i))
    return protid2phyid

def get_seed_or_collateral_tree(phyprot, phylome, p):
    """Return tree object for seed or collateral tree (if not seed tree) 
    for the closest seedid.
    """
    bestSeedid = ""
    tdata = p.get_tree(phyprot, phylome)
    #if no tree for protid as seed look for collateral trees
    if not tdata: 
        #get seedids of collateral trees
        seedids = filter(lambda x: x[2]==phylome, p.get_collateral_seeds(phyprot))
        bestSeedid = None
        smallerstDistance = 999999
        bestTdata = None
        for seedid, seedspcode, phyid in seedids:
            if not seedid.startswith("Phy"):
                seedid = "Phy%s_%s"%(seedid, seedspcode)
            #tdata = p.get_best_tree(seedid, phylome)
            tdata = p.get_tree(seedid, phylome)
            if not tdata: 
                continue 
            bestMethod = sorted(tdata, key=lambda x: tdata[x]['lk'], reverse=True)[0]
            t = ete2.PhyloTree(tdata[bestMethod]['tree'], \
                               sp_naming_function=_get_spcode)
            d = t.get_distance(seedid, phyprot)
            if d < smallerstDistance: 
                smallerstDistance = d
                bestTdata = tdata
                bestSeedid = seedid
        tdata = bestTdata
  
    if not tdata: 
        return None, bestSeedid
    
    #get method giving best lk tree
    bestMethod = sorted(tdata, key=lambda x: tdata[x]['lk'], reverse=True)[0]
  
    #get tree reconstructed with that method
    bestTree = ete2.PhyloTree(tdata[bestMethod]['tree'], \
                              sp_naming_function=_get_spcode)
    
    return bestTree, bestSeedid

def get_orthologs(phyprot, phyid, p, code2score={}):
    """Return orthologs for given protid"""
    orthologs = set((phyprot, ))
    #get tree
    t, bestSeedid = get_seed_or_collateral_tree(phyprot, phyid, p)
    if not t:
        info = "[WARNING] No tree for %s in phylome %s has been found.\n"
        sys.stderr.write(info%(phyprot, phyid))
        return orthologs, code2score, bestSeedid
        
    #set species naming function
    #t.set_species_naming_function(_get_spcode) ###set retrieving species informations!
    #root tree
    if phyid in ROOTED_PHYLOMES:
        outgroup = t.get_farthest_oldest_leaf( ROOTED_PHYLOMES[phyid] )
        t.set_outgroup( outgroup )
    else:
        t.set_outgroup( t.get_midpoint_outgroup() )
    
    #get phyprot node
    l=t.get_leaves_by_name(phyprot)[0]
  
    #get orthologs
    w=0
    for s in filter( lambda e: e.etype=='S', l.get_my_evol_events() ):
        w+=1
        if phyprot in s.in_seqs:
            seqs=s.out_seqs
        else:
            seqs=s.in_seqs
        for o in seqs: 
            orthologs.add( o )
            if _get_spcode(o) not in code2score:
                code2score[_get_spcode(o)] = [ w, ]
            code2score[_get_spcode(o)].append( w )
  
    return orthologs, code2score, bestSeedid

def get_species_in_phylome(phyid, p):
    """Return code2species"""
    code2species={}
    for taxid, code, version, name, source, date in p.get_proteomes_in_phylome(phyid):
        code2species[code]=(taxid, name)
    return code2species

def load_annotation(fn):
    """Return annotation dictionary."""
    prot2ann = {}
    for i, l in enumerate(open(fn), 1):
        ldata = l[:-1].split('\t')
        if l.startswith("#") or len(ldata)<2:
            continue
        prot2ann[ldata[0]] = "\t".join(ldata[1:])
    return prot2ann
  
def profile(handle, out, phyid, protid2phylmeFn, spCode, speciesInRows, \
            annotationFn, verbose):
    """Generate orthologous gene profile."""
    #get phylomeDB connection
    p = _getConnection()
    
    #get protids
    protids=[] #,protid2pfam,protid2change=get_protids( handle,foldChange,foldChangeColumn )
    genes=[]
    for r in SeqIO.parse(handle, 'fasta'): 
        gene = protid = r.id.split('|')[0] #orf19.3038|TPS2
        protids.append(protid)
        #get gene name if present
        if len(r.id.split('|'))>1:
            gene=r.id.split('|')[1].split('_')[0]
        genes.append(gene)
  
    #load annotation
    prot2ann = {}
    if annotationFn:
        prot2ann = load_annotation(annotationFn)       
  
    #get species info
    code2name = get_species_in_phylome(phyid, p)
        
    #define empty profiles
    code2profile={}
    code2score={}
    for code in code2name: 
        code2profile[code]=[0 for i in range(len(protids))]
        code2score[code]=[]
    
    #get phylomedb ids
    protid2seedid = {}
    k = 0
    protid2phyid = get_protid2phyid(protid2phylmeFn, protids, spCode, p, phyid)
    for i, protid in enumerate(protids, 1):
        sys.stderr.write(" %s / %s %s   \r"%(i, len(protids), protid))
        if protid not in protid2phyid: 
            continue
        phyprot = protid2phyid[protid]
        orthologs, code2score, seedid = get_orthologs(phyprot, phyid, p, code2score)
        protid2seedid[protid] = seedid #s; print protid, seedids
        #fill profiles
        for o in orthologs: 
            code2profile[_get_spcode(o)][i-1]+=1
        if len(orthologs)>1:
            k += 1
        elif verbose:
            sys.stderr.write("[WARNING] Only %s orthologs found for %s (%s)!\n"%(len(orthologs), protid, phyprot))
    #write info
    sys.stderr.write("%s proteins; %s with orthologs\n"%(i, k))
        
    ###print summary
    #header
    if not speciesInRows:
        info='#Protid\tGene\tSeedID'
        for code in sorted(code2name, key=lambda x: np.mean(code2score[x])):
            nameShort='%s.%s' % (code2name[code][1][0], code2name[code][1].split()[1])
            info+='\t%s' % nameShort
        info+='\tAnnotation'
        #data
        for j in range(len(protids)):
            protid=protids[j]
            gene=genes[j]
            if gene==protid:
                gene=''
            seedid = ""
            if protid in protid2seedid and protid2seedid[protid]:
                seedid = protid2seedid[protid]
            info+='\n%s\t%s\t%s' % (protid, gene, seedid)
            for code in sorted(code2name,key=lambda x: np.mean(code2score[x])):  
                info+='\t%s' % code2profile[code][j]
      
            if protid in prot2ann:
                info+="\t%s" %  prot2ann[protid] 
    else:
        info='#Species'
        protidLine='#Protid'
        annLine='#Annotation'
        for protid,gene in zip(protids,genes): 
            info+='\t%s' % gene
            if gene!=protid:
                protidLine+='\t%s' % protid
            else:
                protidLine+='\t'
            seedid = ""
            if protid in protid2seedid and protid2seedid[protid]:
                seedid = protid2seedid[protid]
            protidLine += "\t%s" %seedid
            if protid in prot2ann:
                annLine+="\t%s" % prot2ann[protid]
            else:
                annLine+='\t'
        j=0
        for code in sorted(code2name, key=lambda x: np.mean(code2score[x])):  
            info+='\n%s.%s' % (code2name[code][1][0], code2name[code][1].split()[1])
            for orthologNo in code2profile[code]:
                info+='\t%s' % orthologNo
    
        info += protidLine + annLine
    out.write(info)
  
def main():
    import argparse
    usage   = "%(prog)s -v -i CANAL.faa -p63 -s CANAL" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--input",  default=sys.stdin, type=file,
                        help="input stream [stdin]")
    parser.add_argument("-o", "--output",  default=sys.stdout, type=argparse.FileType('w'),
                        help="output stream [stdout]")
    parser.add_argument("-p", "--phyid",  required=True, type=int,
                        help="phylome ID")       
    parser.add_argument("-s", "--spCode", required=True,
                        help="species code")
    parser.add_argument("--speciesInRows", action="store_true", default=False,
                        help="report species in rows [proteins in rows]") 
    parser.add_argument("-a", "--annotation", default=None, type=str,
                        help="annotation file")                 
    parser.add_argument("-c","--protid2phylme", default=None, type=str,
                        help="protid2phylme mappings")
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    #handle .gz
    if o.input.name.endswith('.gz'):
        o.input = gzip.open(o.input.name)
        
    profile(o.input, o.output, o.phyid, o.protid2phylme, o.spCode, o.speciesInRows,\
            o.annotation, o.verbose)

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
