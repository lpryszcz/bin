#!/usr/bin/env python2
"""


USAGE: 

"""

import os, sys
from datetime import datetime
from optparse import OptionParser
from commands import getoutput
from genome_annotation import load_gtf, load_gff

def main():
  """
  """
  usage = "%prog [options] file1 [file2 ... fileN]"
  parser = OptionParser( usage ) 
  
  parser.add_option("-a", dest="annotation", default='', type=str,
                    help="GTF annotation file [%default]" )
  parser.add_option("-f", dest="fdr", default=1e-04, type=float,
                    help="false dicovery rate [%default]" )
  
  ( o, fnames )   = parser.parse_args()
  sys.stderr.write( "Options: %s\nFiles to be processed: %s\n" % ( o, fnames ) )
  
  if not fnames:
    sys.exit( "Speficy at least one input file" )

  prot2ann,prot2locus={},{}
  if o.annotation.endswith(".gff"):
    gene2position,contig2gene=load_gff( o.annotation )
  elif o.annotation:
    gene2position,contig2gene=load_gtf( o.annotation )    
  #process files
  i=0
  samples=[]
  gene2fc={}
  gene2reads={}
  sys.stderr.write( "Processing files...\n" )
  de2sample={}
  for fn in fnames:
    i += 1
    sys.stderr.write( "%s\t%s\n" % ( datetime.now(),fn, ) )
    s = fn.split('.')[0]
    samples.append( s )
    
    #load files
    i=0
    for l in open(fn):
      #skip header
      i+=1
      if i==1:
        continue
        
      id,baseMean,baseMeanA,baseMeanB,foldChange,log2FoldChange,pval,padj = l.split('\t')
      baseMeanA,baseMeanB = float(baseMeanA),float(baseMeanB)
      
      if id not in gene2fc:
        gene2fc[id]=[]
        gene2reads[id]=[] #baseMeanA
      
      #add expression info
      gene2fc[id].append( log2FoldChange )
      gene2reads[id].append( (baseMeanA,baseMeanB) )
      
      try:
        padj=float(padj)
      except:
        continue
      
      if padj<o.fdr:
        if id in de2sample:
          de2sample[id].append( s )
        else:
          de2sample[id]=[s]
      
  #print out  
  header = "#gene\tlocus\t#de\tcontrol"
  for s in samples:
    header += "\t%s log2" % s
  header += "\tannotation"
  print header
  for gene in sorted( de2sample.keys() ):
    #annotation
    function=''
    if gene in gene2position:
      contig,coords,strand,function = gene2position[gene]
      locus="%s:%s-%s %s" % ( contig,coords[0][0],coords[-1][-1],strand )
    
    out = "%s\t%s\t%s" % ( gene,locus,len(de2sample[gene]) )#,gene2control[gene],'\t'.join(gene2fc[gene]),ann )
    #
    i=0
    for reads,fc in zip( gene2reads[gene],gene2fc[gene] ):
      a,b = reads        
      if not i:
          out+="\t%s" % a
      i += 1
      out += "\t%s" % fc
    out += "\t%s" % function
    print out 
      
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
