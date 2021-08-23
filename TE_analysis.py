#!/usr/bin/env python2
###
# Do similarity search (BLAST) of species-specific proteins against TE multifasta in TEdir.
###
import commands
import ete2
import os
import sys
from datetime import datetime 
from Bio import SeqIO

def getFastas( inFpath,p ):
  fastaFpath=inFpath+'.fasta'
  if os.path.isfile( fastaFpath ): return fastaFpath
  print '\tDownloading sequences...',
  i=0; fastas=''
  for line in open(inFpath):
    line=line.strip()
    if not line: continue
    prot=line
    try: 
      seq=p.get_seqid_info(prot)['seq']
      fastas+='>%s\n%s\n' % ( prot, seq )
      i+=1  
    except: print "Cannot get: %s" % prot
    
  outFile=open( fastaFpath,'w' ); outFile.write( fastas ); outFile.close()
  print '%s sequences saved.' % i, datetime.now() 
  return fastaFpath

def runBlast( inFpath,blastFpath,EvalueTh=1e-01,overwrite=False ):
  """Perform tBlastn search (protein query against nucleotide db translated in 6 frames).
  Return blast out file path.
  """
  dbName=blastFpath.split('/')[-1]
  outFpath=inFpath[:-6]+'_%s.blastout' % dbName
  
  if os.path.isfile(outFpath) and not overwrite: return outFpath

  if not os.path.isfile( blastFpath+'.nsq' ): 
    print "\tCreating blast database: %s" % blastFpath
    os.system( 'formatdb -i %s -o F -v 16000 -p F' % blastFpath ) #-p F - dna
  
  print "\tBlasting..."
  blast_cmd = "blastall -p tblastn -FF -m8 -e %s -b 1000 -a 3 -i %s -d %s -o %s" % ( EvalueTh,inFpath,blastFpath,outFpath )
  blast_result=commands.getoutput( blast_cmd )
  if blast_result.startswith('[blastall] FATAL ERROR:'): print blast_cmd; print blast_result
  
  return outFpath

def _getFastas(fastaFpath):
  fastas={}
  for r in SeqIO.parse( open(fastaFpath) ): fastas[r.id]=r.seq
  
  return fastas

def parseBlastOut( outFpath,fastaFpath,EvalueTh,covTh ):
  """
  """
  print "\tAnalysing results..."
  EvalueTh=float(EvalueTh);covTh=float(covTh)
  if covTh: fastas=_getFastas(fastaFpath) # dict of seq by PhylomeDB id
  pQueryid=_line=None; outFile=open( outFpath+'.txt','w' )
  for line in open(outFpath):
    line = line.strip()
    if line.startswith("#") or not line: continue
    try:
      (queryid, hit, identity, alength, mismatches, gaps, qstart, qend, sstart, send, evalue, score) = line.split("\t")
      #Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
    except:
      print line
      continue
    
    evalue,identity=float(evalue),float( identity )
    if evalue>EvalueTh: continue
    if covTh:
      qLen=len( fastas[queryid] )
      overlap=100.0*(1+int(qend)-int(qstart))/qLen 
      if overlap<covTh: continue
    
    if not pQueryid: #executed just once
      _line='%s\t%s' % (queryid,hit)
      pQueryid=queryid
    elif pQueryid!=queryid: 
      outFile.write( _line+'\n' )
      pQueryid=queryid
      _line='%s\t%s' % (queryid,hit)
    else: _line+='\t%s' % hit

  if _line: outFile.write( _line+'\n' ) # save last entry
  outFile.close()

if __name__=='__main__': 
  fastaFpath,dbsFpath,EvalueTh,covTh=sys.argv[1:]
  outFpath=runBlast( fastaFpath,dbsFpath )
  parseBlastOut( outFpath,fastaFpath,EvalueTh,covTh )

