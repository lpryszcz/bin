#!/usr/bin/env python2
#Fetch fragmented assemblies from NCBI
#USAGE: get_fragmented_assemblies.py [taxid]

import os, re, sys, unicodedata
from Bio import Entrez, SeqIO

Entrez.email='lpryszcz@crg.es'
retmax = batchSize = 10000
#rettype='fasta'

def text2fname(value):
    """Convert str to safe fname"""
    value = unicodedata.normalize('NFKD', unicode(value)).encode('ascii', 'ignore')
    value = unicode(re.sub('[^\w\s-]', '', value).strip().lower())
    value = unicode(re.sub('[-\s]+', '_', value).strip().lower())
    return value

#fungi
taxid=4751
if len(sys.argv)>1:
  taxid = int(sys.argv[1])
query="txid%s[Organism:exp] AND 1:100000[Contig N50] NOT 100000:100000000[Scaffold N50]"%taxid

#get assembly ids
handle = Entrez.esearch(db='assembly', term=query, retmax=1000000)
ids = Entrez.read(handle)['IdList']
print "Fetching %s assemblies [ %s ]..." %(len(ids), query)

k = s = 0 
for i, aid in enumerate(ids, 1):
    #get nuccore link
    response = Entrez.read(Entrez.elink(dbfrom="assembly", db="nuccore", id=aid))
    
    #get assembly stats
    assembly_id = response[0]['LinkSetDb'][-1]['Link'][0]['Id']
    handle = Entrez.efetch(db="nuccore", id=assembly_id, rettype="gb")
    r = SeqIO.parse(handle, 'gb').next()
    species = "%s %s"%(r.annotations['organism'], r.name)
    baseName = text2fname(species)
    #write assembly info
    info = open(baseName+'.info', 'w'); info.write(str(r.annotations)); info.close()
    sys.stderr.write(" %s / %s  %s %s %s           \r"%(i, len(ids), species, k, s))
    
    #get contigs
    wgs_id = r.name #r.annotations['wgs'][0]
    #curl http://www.ncbi.nlm.nih.gov/Traces/wgs/?download=ASRE01.1.fsa_nt.gz | zcat
    if not os.path.isfile("%s.contigs.fa"%baseName):
        www = "http://www.ncbi.nlm.nih.gov/Traces/wgs/?download=%s.1.fsa_nt.gz" % wgs_id[:6]
        cmd = "curl -s %s | zcat > %s.contigs.fa" % (www, baseName)
        #print cmd
        os.system(cmd)
        k += 1

    #differentiate between contigs and scaffolds, contigs if only 1 Linkdb
    if len(response[0]['LinkSetDb'])<3 or os.path.isfile("%s.scaffolds.fa"%baseName):
        continue
        
    contigs = [d['Id'] for d in response[0]['LinkSetDb'][0]['Link']]
    #print len(contigs)
    #post NCBI query
    search_handle     = Entrez.epost(db="nucleotide", id=",".join(contigs))
    try:
        search_results    = Entrez.read(search_handle)
    except Exception, e:
        print e
        continue
    webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"]     
    out = open(baseName+'.scaffolds.fa', 'w')
    #fecth all results in batch of batchSize entries at once
    #print "  Fetching %s scaffolds for %s..."%(len(contigs), species)
    for start in range(0, len(contigs), batchSize):
        handle = Entrez.efetch(db="nucleotide", rettype='fasta', retstart=start, 
                               retmax=batchSize, \
                               webenv=webenv, query_key=query_key)
        #report FASTA to output with ignoring empty entries
        out.write("".join(r.format('fasta') for r in SeqIO.parse(handle, 'fasta') if len(r)))
    out.close()
    s += 1
sys.stderr.write("\n%s %s\n"%(k, s))

