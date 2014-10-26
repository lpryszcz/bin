#!/usr/bin/env python
desc="""Retrieve and concatenate transcripts counts from flux capacitor GTF files. 
"""
epilog="""Author:
l.p.pryszcz@gmail.com

12/02/2014 Barcelona
"""

import argparse, gzip, sys
from datetime import datetime

def unload_comment(comment):
	"""Report GTF comment as dict"""
	k2v = {}
	for entry in comment.split('; '):
		k, v = entry.split()
		v = v.strip('"')
		k2v[k] = v
	return k2v

def gtf2counts(i, fn, transcript2counts, geneNames=0):
	""" """
	for line in gzip.open(fn):
		chrname, ctype, feature, s, e, score, strand, space2, comment = line[:-1].split('\t')
		if feature!="transcript":
			continue
		tid = reads = 0
		k2v = unload_comment(comment)
                #note, here gene names may have duplicates, so it's not safe!
                if   geneNames and 'gene_id' in k2v:
                        tid = k2v['gene_id']
		elif 'transcript_id' in k2v:
			tid = k2v['transcript_id']
		if 'reads' in k2v:
			reads = float(k2v['reads'])
		if not tid:
			sys.stderr.write("Warning: File %s: No transcript_id in line: %s\n"%(fn, line))
			continue
		if tid not in transcript2counts:
			transcript2counts[tid] = []
                #combine multiple counts for single gene
                if i+1==len(transcript2counts[tid]):
                        transcript2counts[tid][i]+=reads
                else:
                        transcript2counts[tid].append(reads)
	return transcript2counts

def fluxgtf2counts(fnames, minreads, geneNames, verbose):
        """
        """
        feature = "transcripts"
        if geneNames:
                feature = "genes"
        #load no. of reads for transcripts
        transcript2counts = {}
        sys.stderr.write("Loading flux counts...\n")
        for i, fn in enumerate(fnames, 1):
                sys.stderr.write(" %s / %s  %s          \r"%(i, len(fnames), fn))
                transcript2counts = gtf2counts(i-1, fn, transcript2counts, geneNames)
        
        #report no. of reads for transcripts
        out = sys.stdout
        header = "\t%s\n" % "\t".join(fn.split('.')[0] for fn in fnames)
        out.write(header)
        sys.stderr.write("Storing read counts for %s %s...\n"%(len(transcript2counts), feature))
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
                out.write("%s\t%s\n"%(tid,"\t".join(str(int(round(c))) for c in counts)))

        sys.stderr.write(" %s %s reported\n"%(k, feature))

def main():

    usage   = "%(prog)s [options] -i *.flux.gtf.gz > tcounts.txt" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog )

    parser.add_argument("-v", "--verbose", default=False, action="store_true")
    parser.add_argument(      "--version", default="0.1")
    parser.add_argument("-i", "--input", nargs="+", 
                       help="flux GTF file(s)")
    parser.add_argument("-l", "--minReads", default=10, type=int, 
                       help="min reads per transcript [%(default)s]")
    parser.add_argument("--geneNames", default=False, action="store_true", 
                       help="report genes counts instead of transcripts")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )

    fluxgtf2counts(o.input, o.minReads, o.geneNames, o.verbose)
        
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!   \n")
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
        
