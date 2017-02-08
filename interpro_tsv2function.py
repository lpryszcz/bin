#!/usr/bin/env python
"""Parse InterProScan tsv output and report predicted function for each gene.

cat tsv | python tsv2function.py
"""

import sys

protid2function = {}
for l in sys.stdin:
	ldata  = l.split('\t')
	if len(ldata)<12 or not ldata[12]:
		continue
	protid, function = ldata[0], ldata[12]
	if function.startswith("NULL"):
		continue
	#add prot to dict
	if protid not in protid2function:
		protid2function[protid] = set()
	#add function
	protid2function[protid].add(function)

#decide if split by pipe
pipeSplit = False
pipes_in_id = len(filter(lambda protid: "|" in protid, protid2function.keys()))
#split protid by pipe if not all protids contain pipe
if len(protid2function)>pipes_in_id:
	pipeSplit = True

for protid, functions in sorted(protid2function.items()):
	if pipeSplit:
		for p in protid.split('|'):
			print "%s\t%s" % (p, "; ".join(functions))
	else:
		print "%s\t%s" % (protid, "; ".join(functions))
	
sys.stderr.write("#%s unique proteins annotated with %s functions\n" % (len(protid2function),sum(len(x) for x in protid2function.itervalues())))



