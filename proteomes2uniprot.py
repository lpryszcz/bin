#!/usr/bin/env python
desc="""Crosslink given proteomes and uniprot.
Using md5 hashing to detect the same proteins.
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 6/11/2013
"""

import argparse, os, sys
from Bio      import SeqIO, SwissProt
from datetime import datetime
import hashlib
import gzip

def load_hash2protid(inputs, tposition, tdelimiter, verbose):
    """Return dictionary of taxids and associated md5 hash."""
    hash2protid = {}
    for i, f in enumerate(inputs, 1):
        #get gzip handle
        if f.name.endswith(".gz"):
            f = gzip.open(f.name)
        #get taxid fasta/10.4932.faa.gz -> 4932
        taxid = int(os.path.basename(f.name).split(tdelimiter)[tposition])
        if verbose:
            sys.stderr.write(' %s / %s %s %s        \r' % (i, len(inputs), f.name, taxid))
        #process entries
        hash2protid[taxid] = {}
        for r in SeqIO.parse(f, 'fasta'):
            md5 = hashlib.md5(str(r.seq)).hexdigest()
            if md5 not in hash2protid[taxid]:
                hash2protid[taxid][md5] = [r.id]
            else:
                hash2protid[taxid][md5].append(r.id)
    return hash2protid

def save_entry(outdir, files, protid, db, extid, conf=1.0):
    """Save entry to output file."""
    line = "%s\t%s\t%s\t%s\n" % (protid, db, extid, conf)

    #open file if not yet opened
    db=db.replace('/', '_')
    if not db in files:
        files[db] = gzip.open(os.path.join(outdir,"%s.txt.gz" % db), "w")
        
    files[db].write(line)
    return files
    
def uniprot2metaphors(outdir, hash2protid, accessions, verbose=True):
    """Parser for dat formatted dump of uniprot database. 
    It's quick, much faster than xml parsing!
    But from time to time may encounter wrongly formatted DAT file.
    """
    # create working dir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    #parse dump in dat format
    files = {} # will store opened files
    taxid2stats = {} # { taxid: [matches,total] }
    for r in SwissProt.parse(sys.stdin):
        #skip entry if taxid not of interest
        taxid = int(r.taxonomy_id[0])
        if taxid not in hash2protid: 
            continue
        
        #update stats
        if not taxid in taxid2stats:
            taxid2stats[taxid] = [0, 0]
        taxid2stats[taxid][1] += 1
    
        #check if md5 match
        md5 = hashlib.md5(r.sequence).hexdigest()
        if md5 not in hash2protid[taxid]:
            #skip if no match
            continue 
            
        #save uniprot accession
        for protid in hash2protid[taxid][md5]:
            files  = save_entry(outdir, files, protid, "accessions", r.accessions[0])

            #provide only accessions if requested so 
            if accessions:
                continue

            #add xreference information for each external db
            for ex_db_data in r.cross_references:
                extdb, extid = ex_db_data[:2]
                files = save_entry(outdir, files, protid, extdb, extid)

            #save gene name
            if r.gene_name.startswith('Name='):
                extid = r.gene_name[5:].split(';')[0]
                files = save_entry(outdir, files, protid, "GeneName", extid)

            #save description
            if r.description:
                files = save_entry(outdir, files, protid, "Description", r.description)

            #update stats    
            taxid2stats[taxid][0] += 1

    #close opened files
    for f in files: 
        files[f].close()

    #print stats
    print "#taxid\tmapped\ttotal\t%"
    for taxid in sorted(taxid2stats.keys()):
        mapped,total = taxid2stats[taxid]
        print "%s\t%s\t%s\t%.2f" % (taxid, mapped, total, mapped*100.0/total)
    print "%s\tDone." % datetime.now()
  
def main():
    """
    """
    usage = "wget -O- ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz | zcat | %(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)

    parser.add_argument("-v", "--verbose", default=False, action="store_true")
    parser.add_argument("--version", action='version', version='%(prog)s 1.0')
    parser.add_argument("-i", dest="input", nargs="+", type=file, 
                        help="input fasta(s)")
    parser.add_argument("-o", dest="outdir", default="uniprot",
                        help="output directory [%(default)s]")
    parser.add_argument("-d", dest="db", default='sprot', type=str,
                        help="database         [%(default)s]")
    parser.add_argument("-t", "--tposition", default=1, type=int,
                        help="position of taxid in file name ie fasta/10.4932.faa.txt [%(default)s]")
    parser.add_argument("-u", "--tdelimiter", default=".", 
                        help="delimiter to get taxa [%(default)s]")
    parser.add_argument("--accessions", default=0, type=str,
                        help="map only accessions")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stdout.write( "Options: %s\n" % str(o) )
        
    print "[%s] Loading protein information..." % datetime.ctime(datetime.now())
    hash2protid = load_hash2protid(o.input, o.tposition, o.tdelimiter, o.verbose)
    
    print "[%s] Processing uniprot..." % datetime.ctime(datetime.now())
    outdir = os.path.join(o.outdir, o.db)
    uniprot2metaphors(outdir, hash2protid, o.accessions, o.verbose)
  
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
