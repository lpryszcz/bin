#!/usr/bin/env python
"""Contain functions useful for dealing with phylomeDB."""

import sys
import MySQLdb
#from ete2    import PhylomeDB3Connector as PhylomeDBConnector
from phylomeDB2 import PhylomeDBConnector
from MyGraph import MyGraph
try:
    from rooted_phylomes  import ROOTED_PHYLOMES
except:
    sys.stderr.write("Cannot import ROOTED_PHYLOMES\n")
    

#20131107
PUBLIC_PHYLOMES = [1, 3, 4, 5, 7, 8, 10, 16, 18, 19, 20, 21, 22, 23, 24, 26, 27,
  28, 29, 84, 86, 95, 96, 138, 139, 140, 141, 142, 144, 145, 146, 150, 153, 174,
  179, 183, 184, 191, 137, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
  112, 113, 114, 115, 125, 126, 127, 128, 129, 130, 131, 132, 133, 205, 206,
  207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221,
  222, 223]

def _get_phylomedb_cursor():
    """Return MySQLdb cursor to PhylomeDB."""
    cnx  = MySQLdb.connect(host="mysqlsrv-tgabaldon.linux.crg.es", user="phyReader", \
                           passwd="phyd10.-Reader", db='phylomedb')
    cur  = cnx.cursor()
    return cur

def _getConnection(db="phylomedb", host="mysqlsrv-tgabaldon.linux.crg.es"):
    """Return connection to phylomeDB"""
    p = PhylomeDBConnector(host=host, db=db, user="phyReader", passwd="phyd10.-Reader")
    '''p._trees       = "tree"
    #p._phylomes    = "phylome"
    p._phylomes_table = "phylome"
    #p._algs        = "alignment"
    p._phy_content = "phylome_content"'''
    return p

def _get_spcode(protid):
    """Species naming function compatible with phylome_db3"""
    return protid.split('_')[-1]

def get_evolEvents( t,phyid,seedid ):
    """Root tree and return evol_events"""
    #set species naming function
    t.set_species_naming_function(_get_spcode) 
    #get seed leaf
    try:
        seedNode = t.get_leaves_by_name(seedid)[0]
    except:
        sys.stderr.write( "  Warning: get_evolEvents: Cannot get seed node (%s) in: %s\n" % (seedid,t.write()) )
        return None,None
    #ROOT TREE with error avoinding loops - should be sorted out in the future!    
    try:
        if phyid in ROOTED_PHYLOMES:
            t.set_outgroup( seedNode.get_farthest_oldest_node( ROOTED_PHYLOMES[phyid] ) )
        else:
            t.set_outgroup( t.get_midpoint_outgroup() )
    except:
        sys.stderr.write( "  Warning: get_evolEvents: Cannot root tree: %s %s\n" % (seedid,t.write()) )
        return None,None
      
    #GET EVOLEVENTS
    try:
        evolEvents=t.get_descendant_evol_events()
    except:
        return None,None
  
    return evolEvents,seedNode
    
def get_external( p,seedid,dbs ):
    """Return external ids for given phylome id from requested dbs."""
    extids  = []
    seqinfo =  p.get_seqid_info( seedid )
    #return extids if not external in seqinfo
    if 'external' not in seqinfo:
        return extids
    for db in dbs:
        if db in seqinfo['external']:
          for id in sorted( seqinfo['external'][db] ):
            extids.append( id )
        else:
          extids.append('')
    return extids
    
def get_orthology_graph( t,seedname,phyid ):
    """Return graph representing orthology relationships from the tree."""
    #get evol_events
    evolEvents,seednode=get_evolEvents( t,phyid,seedname )
    if not evolEvents:
        return 
    #create orthologs graph
    ortho_graph=MyGraph( t.get_leaf_names() ) #non-directed graph
  
    #populate it with orthologs  
    for e in filter( lambda x: x.etype=='S',evolEvents ): #seednode.get_my_evol_events() ):
        for o1 in e.out_seqs:
            for o2 in e.in_seqs:
                ortho_graph.add_line( o1,o2 )        

    return ortho_graph
  
def get_age_balanced_outgroup(root, species2age):
    """Root using balanced age
    """
    all_seqs = set(root.get_leaf_names())
    outgroup_dist  = 0
    best_balance = max(species2age.values())
    outgroup_node  = root
    outgroup_size = 0
    sp_fn = root.iter_leaves().next()._speciesFunction
    for leaf in root.iter_descendants():
        leaf_seqs = set(leaf.get_leaf_names())
        size = len(leaf_seqs)
        leaf_species =[sp_fn(s) for s in leaf_seqs]
        out_species = [sp_fn(s) for s in all_seqs - leaf_seqs]
        leaf_age_min = min([species2age[sp] for sp in leaf_species])
        out_age_min = min([species2age[sp] for sp in out_species])
        leaf_age_max = max([species2age[sp] for sp in leaf_species])
        out_age_max = max([species2age[sp] for sp in out_species])
        leaf_age = leaf_age_max - leaf_age_min
        out_age = out_age_max - out_age_min

        age_inbalance = abs(out_age - leaf_age)
        
        # DEBUG ONLY
        # leaf.add_features(age_inbalance = age_inbalance, age=leaf_age)

        update = False
        if age_inbalance < best_balance:
            update = True
        elif age_inbalance == best_balance:
            if size > outgroup_size:
                update = True
            elif size == outgroup_size:
                dist = root.get_distance(leaf)
                outgroup_dist = root.get_distance(outgroup_node)
                if dist > outgroup_dist:
                    update = True

        if update:
            best_balance = age_inbalance
            outgroup_node = leaf
            outgroup_size = size
            
    return outgroup_node
