#!/usr/bin/env python
"""Parse InterProScan tsv output and report GO for each gene.

cat tsv | python tsv2go.py
"""

import sys
import re
gopat = re.compile("GO:\d+")

def extract_go_from_str(gos):
	"""Return GO ids separated by '|'"""
	golist = gopat.findall(gos)
	if golist:
		return "|".join(golist)	
	return ""

protid2go = {}
for l in sys.stdin:
	ldata  = l.split('\t')
	if len(ldata)<14 or not ldata[13]:
		continue
	protid, gos = ldata[0], ldata[13]
	if gos.startswith("NULL"):
		continue
	#extract go terms from str (iprscan4.8)
	if not gos.startswith("GO:"):
		gos = extract_go_from_str(gos)
	if not gos:
		continue
	#add prot to dict
	if protid not in protid2go:
		protid2go[protid] = set()
	#add go
	for go in gos.split('|'):
		protid2go[protid].add(go)
		
for protid, gos in sorted(protid2go.items()):
	print "%s\t%s" % (protid, " ".join(gos))
	
sys.stderr.write("#%s proteins annotated with %s GO terms\n" % (len(protid2go),sum(len(x) for x in protid2go.itervalues())))
