#!/usr/bin/env python
desc="""Fetch entires from NCBI.

Author:
l.p.pryszcz@gmail.com

Barcelona, 15/05/2012
"""

import os, sys
from optparse import OptionParser
from datetime import datetime
from Bio      import Entrez

def ncbi_fetch( queries,taxids,ignored_taxids,db,retmode,batchSize,retmax ):
    """Assign rank to taxid.
    """
    query = '' 
    for q in queries:
        query += "%s OR " % q
    if queries:
        query = "( " + query[:-3] + "AND "
    
    #add required taxids #(txid33208[Organism] OR txid4932[Organism]) NOT (txid7742[Organism] OR txid9606[Organism])
    if taxids:
        query += '('
        for taxid in taxids:
            query+=' txid%s[organism] OR' % taxid
        query = query[:-2] + ')'
  
    #add ignored taxids
    if ignored_taxids:
        query += ' NOT ('
        for taxid in ignored_taxids:
           query+=' txid%s[organism] OR' % taxid
        query = query[:-2] + ')'
  
    #print query
    sys.stderr.write( "Query: %s\n" % query )
  
    #get list of entries for given query
    handle = Entrez.esearch( db=db,term=query,retmax=retmax )
    giList = Entrez.read(handle)['IdList']
  
    #print info about number of proteins
    sys.stderr.write( "Downloading %s entries from NCBI %s database in batches of %s entries...\n" % ( len(giList),db,batchSize ) )
  
    #post NCBI query
    search_handle     = Entrez.epost( db, id=",".join( giList ) )
    search_results    = Entrez.read( search_handle )
    webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"] 
  
    #fecth all results in batch of batchSize entries at once
    for start in range( 0,len(giList),batchSize ):
        #print info
        tnow = datetime.now()
        sys.stderr.write( "\t%s\t%s / %s\n" % ( datetime.ctime(tnow),start,len(giList) ) )
    
        #fetch entries
        error = 1
        #hmmm, this connection catching could be improved
        while error:
            try:
                handle = Entrez.efetch( db=db,retmode=retmode,rettype="txt",retstart=start,retmax=batchSize,webenv=webenv,query_key=query_key )
                fasta  = handle.read()
                error  = 0
            except:
                error += 1
                print "error %s" % error
    
        #print output to stdout
        sys.stdout.write( fasta )

def main():

    usage  = "usage: %prog [options] [ query1 query2 ... queryN ] > out.fa"
    parser = OptionParser( usage=usage,version="%prog 1.0",description=desc ) #allow_interspersed_args=True
    parser.add_option("-d", dest="db",              default='nuccore',
                      help="NCBI database                      [%default]")        
    parser.add_option("-t", dest="taxids",          default='', type=str,
                      help="taxids of interest separated by comma [%default]" )
    parser.add_option("-n", dest="ignored_taxids",  default='', type=str,
                      help="ignored taxids separated by comma  [%default]" )
    parser.add_option("-r", dest="retmode",         default='fasta', type=str,
                      help="output mode: fasta, gb or xml      [%default]" )
    parser.add_option("-e", dest="email",           default="lpryszcz@crg.es", type=str,
                      help="email address                      [%default]" )
    parser.add_option("-m", dest="retmax",          default=100000000, type=int,
                      help="fetch up to m entries in total     [%default]" )         
    parser.add_option("-b", dest="batchSize",       default=10000, type=int,
                      help="fetch b entries at once            [%default]" ) 
    parser.add_option("-v", dest="verbose",  default=False, action="store_true" )      

    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nArgs: %s\n" % ( o,args ) )

    '''for fn in ( o.afile,o.bfile ):
        if not fn:
            parser.error( "Provide input file!" )
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )
    '''
    if not o.email:
        parser.error( "Email (-e) has to be provided!" )
        
    Entrez.email    = o.email
    taxids, ignored_taxids = [], []
    if o.taxids:
        taxids          = o.taxids.split(',')
    if o.ignored_taxids:
        ignored_taxids  = o.ignored_taxids.split(',')
  
    #get from NCBI
    ncbi_fetch( args,taxids,ignored_taxids,o.db,o.retmode,o.batchSize,o.retmax )
          
        
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
