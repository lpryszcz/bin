#!/usr/bin/env python
"""Retrieve one2one orthologs between species A and B, that have been duplicated only once.
So A1 - B1 and A2 - B2 are orthologous, and there is no more orthologs/paralogs between these species.

Require:
ete2

Author:
l.p.pryszcz+git@gmail.com

Barcelona, 4/05/2012
"""

import os, sys
import numpy as np
from optparse import OptionParser
from datetime import datetime
from phylome  import _getConnection,_get_spcode,get_orthology_graph,get_external

def process_tree( p,phyid,n,species_list,seedid ):
  """Return groups of orthologs fulfilling some criteria"""
  groups = t = None
  #get best tree
  trees_dict = p.get_best_tree( seedid,phyid )#; print trees_dict
  #skip if no tree info for given seed
  if not trees_dict:
    return groups,t
  t = trees_dict['tree']
  #skip if no tree    
  if not t:
    return groups,t
      
  #process orthologs as graph
  ortho_graph=get_orthology_graph( t,seedid,phyid )
  if not ortho_graph:
    return groups,t
  #get one2one groups for given set of species
  groups = ortho_graph.nfamily_one2one( species_list,n,_get_spcode )
  return groups,t

def get_consistency( p,t,phyid,seedid,species_list,n,processed_seedids ):
  """Return consistency of given family across
  collateral trees from given phylome."""
  #current tree is ok
  signals = [1]
  processed = set()
  #get seed sp paralogs
  seedsp  = _get_spcode(seedid)
  protids = filter( lambda x: _get_spcode(x)==seedsp,t.get_leaf_names() )
  #process all protids - with iterative update of protids
  while protids:
    #take first protid
    protid = protids.pop()
    #mark as processed
    processed.add( protid )
    #add to processed seedids, so it's visited only once
    processed_seedids.add( protid )
    #skip current tree
    if protid==seedid:
      continue
    #get list of groups A1-B1 A2-B2 and so on
    groups,t = process_tree( p,phyid,n,species_list,protid )
    #add signal
    if groups: #add 1 if correct groups in collateral tree
      signals.append(1)
    else:      #or 0 if collateral tree is inconsistent with tree
      signals.append(0)
    #skip if no t
    if not t:
      continue
    #update protids list - important, as different gene content trees for diverged some families
    for _protid in filter( lambda x: _get_spcode(x)==seedsp,t.get_leaf_names() ):
      if _protid not in processed:
        protids.append( _protid )

  #calculate cs    
  cs = np.mean(signals)
  return cs,processed_seedids
  
def process_phylome( phyid,n,species_list,dbs,step,verbose ):    
  """If not species_list, all species of given phylome are taken."""
  if verbose:
    sys.stderr.write( "[%s] Connecting to PhylomeDB...\n" % datetime.ctime(datetime.now()) )
  
  p=_getConnection()#; print p.get_phylomes() #get some neccesary info
    
  phylome_seedids=p.get_phylome_seed_ids(phyid)[0] #loading seedids #phylome_seedids=['Phy0039MUB_9999994','Phy0039MUC_9999994','Phy0039MQE_9999994']
  if verbose:
    sys.stderr.write( "[%s] Processing %s seeds from phylome_%s...\n" % ( datetime.ctime(datetime.now()),len(phylome_seedids),phyid ) )

  #print header
  header = "#" # "seedid"
  for i in range(n):
    header += "one2one%s\t" % ( i+1, )
  header += "consistency score\n"
  sys.stdout.write( header )

  #process seedids
  i=pI=skipped=positives=0
  processed = set()
  for seedid in phylome_seedids:
    i += 1
    if seedid in processed:
      skipped += 1
      continue

    #get list of groups A1-B1 A2-B2 and so on
    groups,t = process_tree( p,phyid,n,species_list,seedid )
    #do nothing if no such groups
    if not groups:
      continue
      
    #format output line
    line = "" #"%s" % seedid
    ###here you can add protein id conversion
    for group in groups:
      #update processed
      for protid in group:
        processed.add( protid )
        extids = []
        if dbs:
          extids = get_external( p,protid,dbs )
        line += "|".join( [protid,]+extids ) + ","
      line = line[:-1] + "\t"

    #get consistency across collateral trees
    cs,processed = get_consistency( p,t,phyid,seedid,species_list,n,processed )      
    line += "%.3f\n" % cs
    #write to stdout
    sys.stdout.write( line )
    positives += 1
    
    #print progress
    if i>pI:
      pI+=step
      sys.stderr.write( "  %s / %s\t%s    \r" % ( i,len(phylome_seedids),positives ) )
      
  if verbose:
    sys.stderr.write( "[%s] Processed %s seed proteins (duplicated skipped: %s ). %s homologous groups printed.\n" % ( datetime.ctime(datetime.now()),len(phylome_seedids),skipped,positives ) )

def main():
    usage = "usage: %prog [options]" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True

    parser.add_option("-p", dest="phyid", default=0, type=int,
                      help="phylome id               [%default]")
    parser.add_option("-s", dest="species",  default="", 
                      help="comma separated list of species")
    parser.add_option("-n", dest="familySize", default=2, type=int,
                      help="family size ie. for n=2 for species: A, B, C; it requires exactly: A1-B1-C1,A2-B2-C2 [%default]")    
    parser.add_option("-d", dest="dbs",  default="", 
                      help="comma separated list of databases from which add ids [print only phylome id]")
    parser.add_option("-t", dest="step",     default=100, type=int,
                      help="print info every t trees [%default]")
    parser.add_option("-v", dest="verbose",  default=True, action="store_false")
    
    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % ( str(o), ) )

    species = o.species.split(",")
    dbs     = o.dbs.split(",")
    
    if len(species) < 2:
      parser.error( "Provide at least two species codes!" )
    if not o.phyid:
      parser.error( "Provide phylome id!" )
    if o.familySize < 1:
      parser.error( "Family size has to be > 0!" )
      
    process_phylome( o.phyid,o.familySize,species,dbs,o.step,o.verbose )
                      
if __name__=='__main__': 
  t0=datetime.now()
  main()
  #test()  
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
