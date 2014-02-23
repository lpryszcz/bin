#!/usr/bin/env python
desc="""Fetch entires from NCBI.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona/Mizerow, 23/02/2014
"""

import argparse, os, sys
from datetime import datetime
from Bio      import Entrez

def ncbi_fetch(queries, taxids, ignored_taxids, db, rettype, batchSize, retmax, \
               queryAdd, verbose):
    """Fetch from genbank.
    """
    query = '' 
    if queries:      
        query = "(" + " OR ".join(str(q) for q in queries) + ") AND"
    
    #add required taxids #(txid33208[Organism] OR txid4932[Organism]) NOT (txid7742[Organism] OR txid9606[Organism])
    if taxids:
        query += ' (' + " OR ".join("txid%s[organism]"%t for t in taxids) + ')'
    
    #add ignored taxids
    if ignored_taxids:
        query += ' NOT (' + " OR ".join("txid%s[organism]"%t for t in ignored_taxids) + ')'

    #add last bit of query
    if queryAdd:
        query += queryAdd
  
    #print query
    sys.stderr.write( "Query: %s\n" % query )
  
    #get list of entries for given query
    handle = Entrez.esearch(db=db, term=query, retmax=retmax)
    giList = Entrez.read(handle)['IdList']

    #print info about number of proteins
    info = "Downloading %s entries from NCBI %s database in batches of %s entries...\n"
    sys.stderr.write(info % (len(giList), db, batchSize))
  
    #post NCBI query
    search_handle     = Entrez.epost(db, id=",".join(giList))
    search_results    = Entrez.read(search_handle)
    webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"] 
  
    #fecth all results in batch of batchSize entries at once
    for start in range(0, len(giList), batchSize):
        #print info
        tnow = datetime.now()
        sys.stderr.write("[%s] %s / %s\n" % (datetime.ctime(tnow), start+1, len(giList)))
        #fetch entries
        error = 1
        #hmmm, this connection catching could be improved
        while error:
            try:
                handle = Entrez.efetch(db=db, rettype=rettype, retmode="text", \
                                       retstart=start, retmax=batchSize, \
                                       webenv=webenv, query_key=query_key)
                fasta  = handle.read()
                error  = 0
            except:
                error += 1
                sys.stderr.write(" error %s" % error)
        #print output to stdout
        sys.stdout.write(fasta)

def main():

    usage  = "%(prog)s [options] > out.fa"
    parser = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)

    parser.add_argument("-v", "--verbose",  default=False, action="store_true")      
    parser.add_argument("--version", action="version", default="%(prog)s 0.1b")      
    parser.add_argument("--gi",                   nargs="*", type=int, 
                        help="GI to fetch")
    parser.add_argument("-q", "--query",          default="", 
                        help="add str to NCBI query")
    parser.add_argument("-d", dest="db",          default='nuccore',
                        help="NCBI database       [%(default)s]")        
    parser.add_argument("-t", "--taxids",         nargs="+", type=int,
                        help="taxids of interest")
    parser.add_argument("-n", "--ignored_taxids", default='', type=str,
                        help="taxid(s) to ignore  [%(default)s]" )
    parser.add_argument("-e", "--email",          default="lpryszcz@crg.es", 
                        help="email address       [%(default)s]" )
    parser.add_argument("-r", "--rettype",        default='fasta',
                        choices=['fasta', 'gb', 'xml'],
                        help="output mode:        [%(default)s]" )
    parser.add_argument("-m", "--retmax",         default=100000000, type=int,
                        help="max to fetch        [%(default)s]" )         
    parser.add_argument("-b", "--batchSize",      default=10, type=int,
                        help="fetched at once     [%(default)s]" ) 
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n" % o)
  
    #get from NCBI
    Entrez.email = o.email
    ncbi_fetch(o.gi, o.taxids, o.ignored_taxids, o.db, o.rettype, o.batchSize, \
               o.retmax, o.query, o.verbose)
        
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
