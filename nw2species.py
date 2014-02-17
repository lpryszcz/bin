#!/usr/bin/env python 
"""
"""

import os,sys
from ete2 import PhyloTree
from datetime import datetime

def _get_spcode( leaf ):
  return leaf.split('_')[-1]

def main():
  fn=sys.argv[1]
  nw=open(fn).readline()
  
  species={}
  t=PhyloTree(nw)
  
  #set species naming function
  t.set_species_naming_function(_get_spcode)
  
  for l in t.get_leaves():
    spCode=l.species
    try:    species[spCode]+=1
    except: species[spCode] =1
    
  for spCode in sorted( species, key=lambda x: species[x], reverse=True ):
    print '%s\t%s' % ( spCode,species[spCode] )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
