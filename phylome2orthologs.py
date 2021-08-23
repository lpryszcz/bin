#!/usr/bin/env python2
"""Return orthologous groups for given species from a phylome.

Require rooted_phylomes.py (cgenomics.crg.es:/home/services/web/phylomedb.org/wsgi/rooted_phylomes.py)

Based on:
~/workspace/TreeOfLife/src/orthologous_groups.py

Command:
./phylome2orthologs.py -p 120 -s 9999991,9999992,9999993,9999998,CANPA

Author:
l.p.pryszcz+git@gmail.com

Barcelona, 19/04/2012
"""

import os, sys
from optparse import OptionParser
from copy import copy
from datetime import datetime
from rooted_phylomes import ROOTED_PHYLOMES
from ete2 import PhyloTree, PhylomeDB3Connector
from MyGraph import MyGraph
from phylome import _getConnection, _get_spcode, get_evolEvents

def get_orthologs( t,seedname,phyid,species_list ):
    """
    Takes rooted tree as input.
    Return two list: orthology_list - all2all orthologous sequences retrieved from given tree. 
    The group contain one species (the first possible) from each phylogroup.
    """
    #get evol_events
    evolEvents,seednode=get_evolEvents( t,phyid,seedname )
    if not evolEvents:
        return []
  
    orthologs = [ seedname ]
    cur_species=[_get_spcode(seedname)]
    for e in filter( lambda x: x.etype=='S', seednode.get_my_evol_events() ):
        for o in e.out_seqs: 
            spCode = _get_spcode(o)
            if spCode in species_list: 
                orthologs.append( o )
                if spCode not in cur_species:
                    cur_species.append(spCode)
        if len(cur_species)==len(species_list):
            break 
    return orthologs

def populate_orthology_graph( t,seedname,phyid,species_list ):
    #get seed orthologs
    orthologs=get_orthologs( t,seedname,phyid,species_list )
    if not orthologs:
        return [],None
    #create orthologs graph
    ortho_graph=MyGraph( orthologs ) #non-directed graph
    #populate it with orthologs
    for orth in orthologs:
        if orth==seedname:
            continue
        ortho_graph.add_line( seedname, orth )
        #mark all orthologs of given ortholog in the graph
        for orth_in in get_orthologs( t,orth,phyid,species_list ):
            if orth_in in orthologs:
                ortho_graph.add_line( orth, orth_in )
      
    return orthologs,ortho_graph

def generate_orthogroups( orthologs,species_list,collpase_inparalogs=True ): 
    """
    """ 
    sp2orth={}
    for o in orthologs:
        sp=_get_spcode(o)
        try:
            sp2orth[sp].append(o)
        except:
            sp2orth[sp]=[o]
  
    orthogroups=[]
    for sp in species_list:
        if sp in sp2orth:
            #collapse in-paralogs
            if orthogroups and collpase_inparalogs: 
                for o in sp2orth[sp]:
                    orthogroups[-1].append(o)
            #add another species orthologs
            elif orthogroups: 
                orthogroups_org=copy(orthogroups)
                #multiply o.groups when in-paralogs found
                for i in range(len(sp2orth[sp])-1): 
                    for orthogroup in orthogroups_org:
                        orthogroups.append(copy(orthogroup))#; print orthogroups, plen
                index=0
                for o in sp2orth[sp]:
                    for i in range( len(orthogroups)/len(sp2orth[sp]) ): 
                        try:
                            orthogroups[index].append(o)
                        except:
                            print "Error! index: %s\n%s\n%s" % ( index,orthologs,orthogroups )
                    index+=1#make sure adding only one of given species in-paralogs into orthogroups
            #populate orthogroups with first orthologs
            else: 
                for o in sp2orth[sp]:
                    orthogroups.append([o])

    return orthogroups

def get_orthogroups( t,seedname,phyid,species_list,one2one=False,collpase_inparalogs=True ):
    """Looks for orthogroup that is group of all2all orthologs from given set of species (phylogroups)."""
    complete_orthogroups=[]
    #get orthologs and graph representing relationships of orthologs
    orthologs,ortho_graph=populate_orthology_graph( t,seedname,phyid,species_list )
    if not orthologs: return complete_orthogroups
    #get all possible orthologs between species (one protein from each species)
    orthogroups=generate_orthogroups( orthologs,species_list,collpase_inparalogs )
    
    for orthogroup in orthogroups:
        if one2one and ortho_graph.unique( orthogroup ): 
            complete_orthogroups.append( orthogroup ) #if one2one orthologs approve only for connections between orthogroups
        elif not one2one:
            if collpase_inparalogs and ortho_graph.complete_with_inparalogs( orthogroup,_get_spcode ): 
                complete_orthogroups.append( orthogroup )
            elif not collpase_inparalogs and ortho_graph.complete( orthogroup ): 
                complete_orthogroups.append( orthogroup ) #otherwise check only for complete penetration of orthogroup
                
    return complete_orthogroups

def process_phylome( phyid,species_list=None,one2one=True,collpase_inparalogs=False,missingSpeciesTh=0.10,step=100 ):    
    """If not species_list, all species of given phylome are taken.
    """
    print "Generating orthogroups..."
    all_orthogroups=[]
    p=_getConnection()#; print p.get_phylomes() #get some neccesary info
    if not species_list:
        species_list=[]
        proteomes_in_phylome=p.get_proteomes_in_phylome(phyid)['proteomes']
        for proteomeID in proteomes_in_phylome: 
            spCode=proteomeID.split('.')[0]
            species_list.append(spCode)

    print " for %s species: %s" % ( len( species_list ),", ".join( species_list ) ) 
    #make sure seed species if in orthogroups
    seed_sp = p.get_phylome_info(phyid)['seed_proteome'].split('.')[0]
    if not seed_sp in species_list:
        species_list.append( seed_sp )
              
    orthoFpath='phylome%s_orthogroups_%s_%s.txt' % ( phyid,len(species_list),missingSpeciesTh )#; uncommonFpath='phylome%s_uncommon_%s_%s.txt' % ( phyid,len(species_list),missingSpeciesTh )
    if os.path.isfile( orthoFpath ): 
        print " Loading orthologous groups from file: %s" % orthoFpath
        for line in open(orthoFpath):
            line=line.strip()
            all_orthogroups.append( line.split('\t') )
        return all_orthogroups,orthoFpath,species_list
  
    phylome_seedids=p.get_phylome_seed_ids(phyid)[0] #loading seedids #phylome_seedids=['Phy0039MUB_9999994','Phy0039MUC_9999994','Phy0039MQE_9999994']
    outFile=open( orthoFpath,'w' )#; uncommonFile=open( uncommonFpath,'w' )
    trees=not_included=i=pI=low_species_cov=0; pt=datetime.now()
    for seedid in phylome_seedids:
        trees_dict=p.get_best_tree( seedid,phyid )#; print trees_dict
        if not trees_dict:
            continue
        t=trees_dict['tree']
        if not t:
            continue
        trees+=1
        #process orthogroups
        orthogroups=get_orthogroups( t,seedid,phyid,species_list,one2one,collpase_inparalogs )
        for og in orthogroups:
            line=""; _curSpecies=[]
            for o in og: 
                spCode=_get_spcode(o)
                if spCode not in _curSpecies:
                    _curSpecies.append(spCode)
                    line+="%s\t" % o                    
            species_coverage=len(_curSpecies)*1.0/len(species_list)
            line=line[:-1]+'\n'
            #if not enough species in orthogroup 
            if species_coverage<1-missingSpeciesTh:  
                #uncommonFile.write( line) # save in uncommon
                low_species_cov+=1
            else: 
                outFile.write( line ) #save in orthogroup file
                all_orthogroups.append(og) #and add orthogroup to list
                i+=1

        if trees>pI:
            pI+=step
            sys.stdout.write( "   %s %s %s\t%s\r" % ( trees,i,seedid,datetime.now()-pt ) )
            pt=datetime.now()
    print
    print " Processed %s trees (skipped: %s ) for %s seeds. %s one2one orthologous groups and %s with species coverage < %s." % ( trees,not_included,len(phylome_seedids),i,low_species_cov,1-missingSpeciesTh )
    outFile.close()
    return all_orthogroups,orthoFpath,species_list

def main():
    usage = "usage: %prog [options]" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True

    parser.add_option("-p", dest="phyid", default=0, type=int,
                      help="define phylome id                 [%default]")
    parser.add_option("-s", dest="species",  default="", 
                      help="comma separated list of species   [all]")
    parser.add_option("-f", dest="fraction", default=1.0, type=float,
                      help="fraction of species in orthogroup [%default]")
    parser.add_option("-i", dest="one2one",  default=True,  action="store_false",
                      help="allow in-paralogs")
    parser.add_option("-j", dest="collapse", default=False, action="store_true",
                      help="report only one in-paralog orthogroup")
    parser.add_option("-t", dest="step",     default=100, type=int,
                      help="print info every t trees          [%default]")
    parser.add_option("-v", dest="verbose",  default=True, action="store_false")
    
    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "%s\n" % ( str(o), ) )

    missingSpeciesTh = 1.0 - o.fraction
    species = o.species.split(",")
    if len(species)<2:
        species = []
    #print species
    process_phylome( o.phyid,species,o.one2one,o.collapse,missingSpeciesTh,o.step )
  
if __name__=='__main__': 
    t0=datetime.now()
    main()
    dt=datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
