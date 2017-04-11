#!/usr/bin/env python
"""
Download genomes from NCBI for given taxid.

USAGE: get_complete_genomes.NCBI.py taxid nuccore
"""

import os, sys
from datetime import datetime
from Bio      import Entrez, SeqIO
Entrez.email = 'lpryszcz@crg.eu'

def split_seq( seq,length=70 ):
  """Return list of sequence fragnments having given length.
  """
  return [seq[i:i + length] for i in range(0, len(seq), length)]
 
def gb2cds( handle,outDir='out' ):
  """
  """
  #get genome sequence as GenBank SeqIO object
  for gb in SeqIO.parse( handle,"gb" ): #Entrez.efetch( db=db,id=gi,rettype="gb" )
    ###get sp name
    acc   = gb.id
    spName= gb.description.split(',')[0]
    spName= spName.strip().replace(' ','_').replace('/','_').replace(':','_')
    
    #check if outfile exists
    outFn = '%s/%s.%s.chromosome.fna' % (outDir,spName,acc)
    #if os.path.isfile( outFn ):
    #  continue
    
    ###save genome fasta
    seq   = '\n'.join( split_seq( str(gb.seq) ) )
    out   = open( outFn,'w' )
    out.write( '>%s\n%s\n' % ( acc,seq ) )
    out.close()
    
    ###save genes, CDS and translations
    #genesOut = open( '%s/%s.%s.genes.fna' % (outDir,spName,acc),'w' )
    cdsOut   = open( '%s/%s.%s.cds.fna'   % (outDir,spName,acc),'w' )
    aaOut    = open( '%s/%s.%s.aa.faa'    % (outDir,spName,acc),'w' )
    for feature in gb.features:
      if feature.type == 'CDS':
        #header
        header = '>'
        if 'db_xref' in feature.qualifiers:
          extID = feature.qualifiers['db_xref'][0]
          extID = extID.replace(':','|')
          header += "%s|" % extID
        
        if 'protein_id' in feature.qualifiers:
          header += "%s"    % feature.qualifiers['protein_id'][0]
        header += '|'
        
        if 'gene' in feature.qualifiers:
          header += ' %s '   % feature.qualifiers['gene'][0]
        header += '|'
        
        #add position info
        header += ' %s-%s %s' % ( feature.location.start.position,feature.location.end.position,feature.strand )
        
        #cds
        cdsSeq  = feature.extract( gb.seq )
        cds     = '\n'.join( split_seq( str( cdsSeq ) ) )
        cdsOut.write( '%s\n%s\n' % ( header,cds) )
        
        #protein
        if 'translation' in feature.qualifiers:
          aa  = '\n'.join( split_seq( feature.qualifiers['translation'][0] ) )
        else:
          aa  = '\n'.join( split_seq( str( cdsSeq.translate() ) ) )
        aaOut.write( '%s\n%s\n'  % ( header,aa ) )
        
    #close out files
    cdsOut.close()
    aaOut.close()

def get_cds( outDir = 'cds' ):
  """Get protein and cds sequence for every chromosome
  """
  #create outdir
  if not os.path.isdir( outDir ):
    os.makedirs( outDir )
  
  print "Retrieving CDS & proteins from GenBank files..."  
  fnames = filter( lambda x: x.endswith('.gb'), os.walk('.').next()[2] )
  fnames.sort()
  for fn in fnames:
    tNow=str(datetime.now())  
    print "%s\t%s / %s\t%s" % ( tNow.split('.')[0],fnames.index(fn)+1,len(fnames),fn )
    gb2cds( open(fn),outDir )

def get_genomes( taxid,db ):
  """Fetch complete genomes in GenBank format
  """
  q='txid%s[organism] AND biomol_mrna[PROP]' % taxid #'txid%s[organism] AND complete genome[title]' % taxid 
  handle = Entrez.esearch( db=db,term=q,retmax=1000000 )
  record = Entrez.read( handle )
  GIs=record['IdList']
  
  print "Downloading %s complete genomes...\n cmd: %s" % ( len(GIs),q )
  #fecth all genomes
  for gi in GIs:
    #skip if file already exits
    outFn = '%s.gb' % gi
    if os.path.isfile( outFn ):
      continue
    
    #print info  
    tNow=str(datetime.now())  
    print "%s\t%s / %s\t%s" % ( tNow.split('.')[0],GIs.index(gi)+1,len(GIs),gi )
    
    #get GenBank string and write
    error = 1
    while error:
      try:
        gb     = Entrez.efetch( db=db,id=gi,rettype="gb" ).read()
        error  = 0
      except:
        error += 1
        print "\terror %s" % error

    #skip fragmented entries like: http://www.ncbi.nlm.nih.gov/nuccore/AL935263.1
    #if len(gb)<10**5:
    #  continue
    
    #open outfile, write, close
    out   = open( outFn,'w' ); out.write( gb ); out.close()
    
    #remove wrong entries
    #try:
    #  r = SeqIO.parse( open(outFn),'gb' ).next()
    #except:
    #  print " Wrong genbank file: %s. Removed!" % outFn 
    #  os.unlink( outFn )
    
def main(groupTaxid, db="nuccore"):
  """
  """
  ###get taxIDs of strains with complete genome, groupTaxid=2 for Bacteria
  get_genomes( groupTaxid,db )
  
  ###get CDSs and proteins
  get_cds()

if __name__=='__main__': 
  t0=datetime.now()
  groupTaxid = int(sys.argv[1])
  db=sys.argv[2]
  main(groupTaxid, db)
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
