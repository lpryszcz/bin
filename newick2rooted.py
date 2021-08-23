#!/usr/bin/env python2

import sys
import ete2
rooters = False
try:
	from rooted_phylomes import ROOTED_PHYLOMES
	rooters = True
except:
	sys.stderr.write("No rooters (rooted_phylomes.py)\n")
try:	
	phyid, seed = int(sys.argv[1]), sys.argv[2]
except:
	sys.stderr.write("No phylome_id and/or seed given\n")
	rooters = False

def _get_spcode(protid):
    """Species naming function compatible with phylome_db3"""
    return protid.split('_')[-1]

for i, nw in enumerate(sys.stdin, 1):
	if not i%100:
		sys.stderr.write(" %i   \r"%i)
	t = ete2.PhyloTree(nw)
	t.set_species_naming_function(_get_spcode) 
	try:
		seedNode = t.get_leaves_by_name(seed)[0]
		seedNode.get_farthest_oldest_node(ROOTED_PHYLOMES[phyid])
		if rooters:
			sys.stderr.write("[WARNING] Cannot root tree %s with rooter\n"%i)
	except:
		t.set_outgroup(t.get_midpoint_outgroup())
	print t.write(format=9)


