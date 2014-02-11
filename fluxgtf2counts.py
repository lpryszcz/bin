#!/usr/bin/env python
"""Retrieve and concatenate transcripts counts from flux capacitor GTF files

USAGE:
python fluxgtf2counts.py *.flux.gtf > transcripts.counts.txt
"""

import gzip, sys
fnames   = sys.argv[1:]
minreads = 10

def unload_comment(comment):
	"""Report GTF comment as dict"""
	k2v = {}
	for entry in comment.split('; '):
		k, v = entry.split()
		v = v.strip('"')
		k2v[k] = v
	return k2v

def fluxgtf2counts(fn, transcript2counts):
	""" """
	for line in gzip.open(fn):
		chrname, ctype, feature, s, e, score, strand, space2, comment = line[:-1].split('\t')
		if feature!="transcript":
			continue
		tid = reads = 0
		k2v = unload_comment(comment)
		if 'transcript_id' in k2v:
			tid = k2v['transcript_id']
		if 'reads' in k2v:
			reads = int(round(float(k2v['reads'])))
		if not tid:
			sys.stderr.write("Warning: File %s: No transcript_id in line: %s\n"%(fn, line))
			continue
		if tid not in transcript2counts:
			transcript2counts[tid] = []
		transcript2counts[tid].append(reads)
	return transcript2counts
	
#load no. of reads for transcripts
transcript2counts = {}
sys.stderr.write("Loading flux counts...\n")
for i, fn in enumerate(fnames, 1):
	sys.stderr.write(" %s / %s  %s     \r"%(i, len(fnames), fn))
	transcript2counts = fluxgtf2counts(fn, transcript2counts)
	
#report no. of reads for transcripts
out = sys.stdout
header = "\t%s\n" % "\t".join(fn.split('.')[0] for fn in fnames)
out.write(header)
sys.stderr.write("Storing read counts for %s transcripts...\n"%len(transcript2counts))
k = 0
for tid, counts in transcript2counts.iteritems():
	#check if enough counts
	if len(counts)<len(fnames):
		sys.stderr.write("Warning: Not enough (%s) read counts for %s\n"%(len(counts), tid))
		continue
	#check if enough reads for transcript
	if sum(counts)<minreads:
		continue
	k += 1
	out.write("%s\t%s\n"%(tid,"\t".join(str(c) for c in counts)))
	
sys.stderr.write(" %s transcripts reported\n"%k)
