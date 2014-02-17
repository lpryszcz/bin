#!/usr/bin/env python
"""Annotates genes having mapping to PFAM from given species with GO codes.
Get http://www.geneontology.org/external2go/pfam2go
"""

import os, sys
from datetime import datetime
from optparse import OptionParser
from Bio      import SeqIO
from pfam     import load_pfam2go, load_broad_pfam, load_pfam_tblout

def gene2go( pfam2goFn,fn,tblout,eTh=1e-05,source='broad'):
  """
  """
  pfam2go = load_pfam2go( pfam2goFn )
  
  if   source=="broad":
    prot2pfam = load_broad_pfam( tblout,eTh )
  elif source=="tblout":
    prot2pfam = load_pfam_tblout( tblout,eTh )
  else:
    sys.exit("Specify valid source: tblout or broad!")
  
  out = open( fn+'.go','w')
  for prot in prot2pfam:
    j=0
    line='%s\t' % prot
    for pfam in prot2pfam[prot]:
      pfamShort = pfam.split('.')[0]
      if not pfamShort in pfam2go:
        continue
      for go in pfam2go[pfamShort]:
        line+='%s;' % go
        j+=1
        
    if j:
      out.write( line.strip(';')+"\n" )
  out.close()
  
def main():
  parser = OptionParser() #allow_interspersed_args=True
  
  parser.add_option("-i", dest="ids", default='', 
                    help="ids input file             [%default]")
  parser.add_option("-t", dest="tblout", default='', 
                    help="hmmer tblout file          [%default]")
  parser.add_option("-p", dest="pfam2goFn", default='pfam2go', 
                    help="pfam2go association file   [%default]" )
  parser.add_option("-e", dest="eTh", default=1e-05, type=float,
                    help="pfam mapping E-value cut-off [%default]")
  parser.add_option("-s", dest="source", default='tblout', 
                    help="mapping source: tblout or broad [%default]" )
  parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False)                    
  
  ( o, args )   = parser.parse_args()
  if o.verbose:
    sys.stderr.write( "Options: %s\nArgs: %s\n" % ( o,args ) )

  if not o.tblout:
    o.tblout = o.ids+".tblout"
  for fn in (o.pfam2goFn,o.ids,o.tblout):
    if not fn:
      parser.error("No all parameters specified!")
    if not os.path.isfile(fn):
      parser.error("No such file: %s")
  
  if o.verbose:
    sys.stderr.write( "Generating association file > %s\n" % (o.ids+".go",) )
  gene2go( o.pfam2goFn,o.ids,o.tblout,o.eTh,o.source )
  
  '''if o.verbose:
    sys.stderr.write( "Generating population file > %s\n" % (o.ids+".ids",) )
  ids = set()
  for r in SeqIO.parse(o.ids,"fasta"):
    prot = r.id
    #correct phylomedb isoforms Phy003WV3G_9999994_2 > Phy003WV3G_9999994
    if prot.startswith("Phy") and prot.count("_")>1:
      prot = "_".join(prot.split("_")[:2])
    ids.add( prot )
    
  out = open(o.ids+".ids","w"); out.write( "\n".join( ids ) ); out.close()'''
    
if __name__=='__main__': 
  t0  = datetime.now()
  main()
  dt  = datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
