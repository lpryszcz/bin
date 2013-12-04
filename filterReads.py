#!/usr/bin/env python
"""
Filter QSEQ (GERALD)/FastQ reads. Store output as FastQ.
Reads are clipped at first undetermined base (. in qseq or N in fastq)
and at first base having qual below -q.
Reads (and their pairs if -p) not passing filtering are discarded. 
Orphaned reads may be store optionally (-u).

USAGE:
filterReads.py -l31 -q10 -o outDir -p -u -v _archives/wt_10b_read*

AUTHOR:
Leszek Pryszcz
lpryszcz@crg.es

Version 0.11

Ideas:
-add automatic qual_offset recognition

Fixes:
-0.1
--wrong sep in _clipSeq

-0.11
--output always PHRED+33 quals (Sanger, CASAVA1.8+)
--include reads with '.' bases -> 'N'
"""
 
import gzip, os, sys
from datetime import datetime
from Bio import SeqIO
from optparse import OptionParser
import locale
locale.setlocale(locale.LC_ALL, 'en_US.utf8')

def _clipSeq( seq,quals,minLen,sep='.' ):
  """
  """
  if sep in seq:
    pos=seq.index(sep)
    if pos<minLen:
      return
    seq,quals=seq[:pos],quals[:pos]
  return seq
  
def _getFastQ( file,qseq,minLen,qualityTh,ASCII_offset,simpleHeaders,fasta,seqi,pair='' ):
  """Process each line of GERALD output and Return FastQ str. Cut seq @ first '.' position.
  Also, check for quality if qualityTh defined.
  Return None if length of the seq < minLen or an error occured during line processing.
  """
  ##GERALD (QSEQ)
  if qseq:
    qseq=file.next()[:-1]
    
    qseq_element=qseq.split('\t') #SOLEXA 90403 4 1 23 1566 0 1 ACCGCTCTCGTGCTCGTCGCTGCGTTGAGGCTTGCG `aaaaa```aZa^`]a``a``a]a^`a\Y^`^^]V` 1
    if len(qseq_element)!=11 or qseq_element[-1]!='1': 
      return
    
    #formatting
    name      = '@%s:%s:%s:%s:%s#%s/%s' % ( qseq_element[0],qseq_element[2],qseq_element[3],qseq_element[4],qseq_element[5],qseq_element[6],qseq_element[7] )
    seq,quals = qseq_element[8],qseq_element[9]
    
    #clip seq & quals @ . ( unknown base )
    seq       = _clipSeq( seq,quals,minLen,'.' )

  ##FASTQ
  else:
    name  = file.next()[:-1]
    seq   = file.next()[:-1]
    sepli = file.next()[:-1]
    quals = file.next()[:-1]
    
    
    #ENA: @ERR005019.1 HWI-EAS111_5_FC13331_PE_R1:8:1:54:539/1 > @ERR005019.1/1
    #format name - @HWI-ST227:145:C06RAACXX:7:1101:1156:2148 2:Y:0:AACT > @HWI-ST227:145:C06RAACXX:7:1101:1156:2148/2
    if len( name.split() ) > 1:
      name = '%s%s' % ( name.split()[0],pair )
    
    #clip seq & quals @ N ( unknown base )
    seq = _clipSeq( seq,quals,minLen,'N' )
  
  if not seq:
    return
  
  #return PHRED+33 quals (Sanger encoding)
  if ASCII_offset!=33:
    quals=''.join( [ chr(ord(q)-ASCII_offset+33) for q in quals ] )
  
  #cut sequence & quals @ quality
  if qualityTh:
    #pos=0
    #for q in quals:
    for pos in range(len(seq)):
      phredQ=ord(quals[pos])-33 #PHRED+33 encoding
      if phredQ<qualityTh:
        seq=seq[:pos]
        break
    #clip seq and qual
    quals=quals[:len(seq)]
  
  if len(seq)<minLen or not seq:  
    return
  
  #define fastQ line
  if simpleHeaders:
    name="@%s%s" % ( seqi,pair )
  if fasta:
    fastq='>%s\n%s\n' % ( name[1:],seq )
  else:
    fastq='%s\n%s\n+\n%s\n' % ( name,seq,quals )
  
  return fastq

def filterSingle( qseq,fPathF,outFileF,minLen,qualityTh,ASCII_offset,\
                  simpleHeaders,fasta,gzlipEndings=('.gz','.tgz','.tar.gz'),step=10**4 ):
  """Convert GERALD file (fPathF) to FastQ (outFileF) with quality and read length filtering.
  Trim @ first base with quality < qualityTh, distard reads with lenth < minLen. 
  """
  #open input file
  i=pI=filterSkip=0
  if fPathF.endswith( gzlipEndings ):
    fileF=gzip.open( fPathF,'rb' )
  else:
    fileF=open( fPathF,'rb' )
  #parse file
  try:
    while 1:
        i+=1
        rec=_getFastQ( fileF,qseq,minLen,qualityTh,ASCII_offset,simpleHeaders,fasta,i )
        
        if rec: 
          outFileF.write(rec)
        else:   
          filterSkip+=1
          
        if i>pI:
          pI+=step
          sys.stderr.write( "%s processed [%.2f%s passed filtering]\r" % ( locale.format("%d",i,grouping=True),(i-filterSkip)*100.0/i,'%' ) )
          #sys.stderr.write( "%s done. %s successes [%.2f%s]\r" % (i,i-filterSkip,(i-filterSkip)*100.0/i,'%') )
          
  except StopIteration, e: 
    pass
  
  return i,filterSkip,0

def filterPaired( qseq,fPathF,fPathR,outFileF,outFileR,combinedOutFile,\
                  outUnpaired,minLen,qualityTh,ASCII_offset,simpleHeaders,\
                  fasta,gzlipEndings=('.gz',),step=10**4 ):
  """Convert GERALD files (fPathF and fPathR) to FastQ (outFileF,outFileR) with quality and read length filtering.
  Trim @ first base with quality < qualityTh, distard reads with lenth < minLen. 
  Write singletons FastQ to outUnpaired if exists.
  Write filtered pairs to combinedOutFile if exists.
  """
  #open input files
  i=pI=filterSkip=both=fori=revi=0
  '''#dealing with remote files
  elif fPathF.startswith("ftp://") and fPathR.startswith("ftp://"):
    wgetF=subprocess.Popen( [ "wget","-O-",fPathF,"|","zcat" ],stdout=subprocess.PIPE )
    wgetR=subprocess.Popen( [ "wget","-O-",fPathR,"|","zcat" ],stdout=subprocess.PIPE )
    fileF=wgetF.stdout
    fileR=wgetR.stdout
  '''
  if fPathF.endswith(gzlipEndings) and fPathR.endswith(gzlipEndings):
    fileF=gzip.open( fPathF,'rb' )
    fileR=gzip.open( fPathR,'rb' )
  else:
    fileF=open( fPathF,'rb' ) 
    fileR=open( fPathR,'rb' )
  #parse files
  try:
    while 1:
        i+=1
        rec1=_getFastQ( fileF,qseq,minLen,qualityTh,ASCII_offset,simpleHeaders,fasta,i,'/1' )
        rec2=_getFastQ( fileR,qseq,minLen,qualityTh,ASCII_offset,simpleHeaders,fasta,i,'/2' )
        
        #save FastQ
        if rec1 and rec2:           #of given pair is F & R pass filtering
          if outFileF:
            outFileF.write(rec1)
            outFileR.write(rec2)
          #store combinedOut if exists
          if combinedOutFile: 
            combinedOutFile.write(rec1+rec2)
          both += 1
        elif outUnpaired and rec1:  #F as single read if R didn't pass filtering
          fori+=1
          outUnpaired.write( rec1 )
        elif outUnpaired and rec2:  #R as single read if F didn't pass filtering
          revi+=1
          outUnpaired.write( rec2 )
        else: 
          filterSkip+=1             #nothing if both didn't pass filtering 
          
        if i>pI:
          pI+=step
          sys.stderr.write( "%9s processed [%6.2f%s ok] Both: %s Single F/R: %s/%s\r" % ( locale.format("%d",i,grouping=True),(i-filterSkip)*100.0/i,'%',both,fori,revi ) )
          #sys.stderr.write( "%s done. %s successes [%.2f%s]\r" % (i,i-filterSkip,(i-filterSkip)*100.0/i,'%') ) #sys.stdout.write( ' %s %s\r' % (i,filterSkip) )
          
  except StopIteration, e: 
    pass
  sys.stderr.write( "%9s processed [%6.2f%s ok] Both: %s Single F/R: %s/%s\n" % ( locale.format("%d",i,grouping=True),(i-filterSkip)*100.0/i,'%',both,fori,revi ) )
  #close in files
  fileF.close()
  fileR.close()
  if combinedOutFile:
    combinedOutFile.close()
    
  return i,filterSkip,fori+revi

def processReads( fPaths,qseq,outDir,paired,storeUnpaired,minLen,qualityTh,\
                  separateFastQ,combinedFastQ,ASCII_offset_33,replace,\
                  simpleHeaders,fasta,verbose,qseqEnds=('.txt','.gz') ):  
  """Convert GERALD/FastQ files from Fpaths list to FastQ with quality and read length filtering.
  Trim @ first base with quality < qualityTh, distard reads with lenth < minLen. 
  If paired, save filtered read pairs into qXX_1.fastq and qXX_2.fastq and remove reads with missing pairs or save them as singletons into qXX.unpaired.fastq (storeUnpaired).
  If combinedFastQ, save F & R into qXX.combined.fastq.
  Print stats at the end.
  """
  #define quality score type
  if ASCII_offset_33: ASCII_offset=33 #sanger
  else:               ASCII_offset=64 #solexa/illumina
  #info
  if verbose: 
    if qseq:
      print "\nConverting gerald file(s) %s into fastq and filtering..." % fPaths, datetime.now()
    else:
      print "\nFiltering FastQ file(s) %s..." % fPaths, datetime.now()
  #process geralds
  i=filterSkip=single=0
  ###
  #process paired-ends reads
  if paired: 
    #define file names for output
    if fasta:
      outFnameF       =os.path.join( outDir,'q%s_1.fasta'         % qualityTh )
      outFnameR       =os.path.join( outDir,'q%s_2.fasta'         % qualityTh )
      unpairedFname   =os.path.join( outDir,'q%s.unpaired.fasta'  % qualityTh )
      outCombinedFname=os.path.join( outDir,'q%s.combined.fasta'  % qualityTh )      
    else:
      outFnameF       =os.path.join( outDir,'q%s_1.fastq'         % qualityTh )
      outFnameR       =os.path.join( outDir,'q%s_2.fastq'         % qualityTh )
      unpairedFname   =os.path.join( outDir,'q%s.unpaired.fastq'  % qualityTh )
      outCombinedFname=os.path.join( outDir,'q%s.combined.fastq'  % qualityTh )
    #check if out file exists
    if replace: pass
    elif os.path.isfile( outFnameF ) or os.path.isfile( outFnameR ):
      print " At least one of out files exist: %s or %s. Exitting." % ( outFnameF,outFnameR )
      return
      
    #open files for writting
    if separateFastQ:
      outFileF=open( outFnameF,'wb' )
      outFileR=open( outFnameR,'wb' )
    else:
      outFileF=outFileR=False
      
    #open out file for unpaired reads
    if storeUnpaired: 
      outUnpaired=open( unpairedFname,'wb' )
    else:             
      outUnpaired=False
      
    #open out file for combined FastQ
    if combinedFastQ: 
      combinedOutFile =open( outCombinedFname,'wb' )
    else:             
      combinedOutFile =False
      
    #store pairs of filenames
    fnames_pair=[]
    #process all input files
    for fname in fPaths:
      fnames_pair.append(fname)#; print fnames_pair
      if len(fnames_pair)!=2: continue
      #get F and R qseq fnames
      fPathF,fPathR=fnames_pair 
      #proces qseq files: GERALD->FASTA
      p_i,p_filterSkip,p_single = filterPaired( qseq,fPathF,fPathR,outFileF,outFileR,combinedOutFile,outUnpaired,minLen,qualityTh,ASCII_offset,simpleHeaders,fasta )
      #print info
      if verbose: 
        print '',fnames_pair,p_i,p_filterSkip,datetime.now()
      #update read counts
      i           +=p_i
      filterSkip  +=p_filterSkip
      single      +=p_single
      #reset fnames
      fnames_pair=[]
      
    #close outfiles
    if outFileF:
      outFileF.close(); outFileR.close()
    #store optional out files
    if storeUnpaired: 
      outUnpaired.close()
    if combinedFastQ: 
      combinedOutFile.close()
    #print info for run
    print '\tProcessed pairs: %s. Filtered: %s. Reads pairs included: %s [%.2f%s]. Singletons: %s [%.2f%s]' % ( i,filterSkip,(i-filterSkip),(i-filterSkip)*100.0/i,'%',single,single*100.0/i,'%' )
    
  ####
  #for single reads
  else: 
    #define out fname
    if fasta:
      outFnameF=os.path.join( outDir,'q%s.fasta' % ( qualityTh ) )
    else:
      outFnameF=os.path.join( outDir,'q%s.fastq' % ( qualityTh ) )
    #check if out file exists
    if replace: pass
    elif os.path.isfile( outFnameF ):
      print " Out file exists: %s. Exitting." % ( outFnameF, )
      return
    #open files for writting
    outFileF=open( outFnameF,'wb' )
    #process all files as single reads
    for fPathF in fPaths: 
      #proces qseq file: GERALD->FASTA
      p_i,p_filterSkip,p_single = filterSingle( qseq,fPathF,outFileF,minLen,qualityTh,ASCII_offset,simpleHeaders,fasta )
      #print info
      if verbose: print '',fPathF,p_i,p_filterSkip
      #update read counts
      i           +=p_i
      filterSkip  +=p_filterSkip
    #close outfile
    outFileF.close()
    #print info for run
    print '\tProcessed reads: %s. Filtered: %s. Reads included: %s [%.2f%s].' % ( i,filterSkip,(i-filterSkip),(i-filterSkip)*100.0/i,'%' )
  
def main():
  
  parser = OptionParser() #allow_interspersed_args=True
  
  parser.add_option("-o", dest="outDir", help="define where to store output files")
  parser.add_option("-v", action="store_true", dest="verbose", default=False,
                    help="print status messages to stdout [default: %default]")
  parser.add_option("-g", dest="qseq", action="store_true",  default=False,
                    help="input QSEQ [default: FastQ]")
  parser.add_option("-l", dest="minLen", default=31, type=int,
                    help="min read lenght (shorter after quality trimming are removed) [default: %default]" )
  parser.add_option("-q", dest="qualityTh", default=0, type=int,
                    help="read is clipped @ first base having PHRED quality lower than [default: %default]" )
  parser.add_option("-t", dest="ASCII_offset_33", action="store_false", default=True,
                    help="use illumina/solexa quality encoding (ASCII offset of 64) [default: Sanger]")
  parser.add_option("-p", action="store_true", dest="paired", default=False, 
                    help="process as paired-end reads (qXX_1.fastq & qXX_2.fastq) [default: single-end]")
  parser.add_option("-u", action="store_true", dest="unpaired", default=False, 
                    help="store orphaned reads > qXX.unpaired.fastq [default: %default]")                  
  parser.add_option("-r", dest="replace", action="store_true", default=False,
                    help="overwrite output files [default: %default]")
  parser.add_option("-b", dest="separateFastQ", action="store_false", default=True,
                    help="don't store separate fastQ for L & R reads > qXX_1.fastq qXX_2.fastq [default: %default]" )                  
  parser.add_option("-c", dest="combinedFastQ", action="store_true", default=False,
                    help="store combined fastQ for paired reads as well > qXX.combined.fastq [default: %default]" )
  parser.add_option("-H", dest="simpleHeaders", action="store_true", default=False,
                    help="replace headers by int [default: %default]" )
  parser.add_option("--fasta", dest="fasta", action="store_true", default=False,
                    help="report fasta           [default: fastq]" )
  
  ( o, fPaths ) = parser.parse_args()
  if o.verbose: 
    print o, fPaths
  print "###\nYou may consider using newer, slighly faster implementation: filterReads.new.py\n###"
  
  if not o.outDir: 
    sys.exit( "Output dir has to be specified.")
  
  ###check input parameters
  #if any file specified
  if not fPaths: sys.exit( "Error! qseq/fastq files has to be specified as input" )
  #file input are files
  for fpath in fPaths:
    if not os.path.isfile( fpath ): 
      sys.exit( "Error! No such file: %s" % fpath )
  
  #create output directore if not present already
  if not os.path.isdir( o.outDir ): 
    os.makedirs( o.outDir )
  
  #gerald2fastq
  processReads( fPaths,o.qseq,o.outDir,o.paired,o.unpaired,o.minLen,o.qualityTh,\
    o.separateFastQ,o.combinedFastQ,o.ASCII_offset_33,o.replace,o.simpleHeaders,\
                o.fasta,o.verbose )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  print "#Time elapsed: %s" % dt
