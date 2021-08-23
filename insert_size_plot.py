#!/usr/bin/env python2
"""

USAGE: for f in *.bam; do echo `date` $f; samtools view $f -f67 -F12 | head -n100000 | insert_size_plot.py $f.insert_size 10000; done

"""
import os,sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from math import fabs

def histplot( sizes,fname ):
  """
  #mu, sigma = 100, 15
  #x = mu + sigma * np.random.randn(10000)
  """
  mu = np.mean( sizes )
  sd = np.std( sizes )
  median = np.median( sizes )
  
  print '%s\tmean = %.2f +- %.2f\t median = %s ( min%s - %s )' % ( fname,mu,sd,median,min(sizes),max(sizes) )
  
  # limit to only 3* mean, usefull for highly dispersed data
  sizes = filter(lambda x: x<3*mu, sizes)
  # the histogram of the data
  n, bins, patches = plt.hist(sizes, bins=max(sizes)/5, normed=1, facecolor='g', alpha=0.75)

  plt.xlabel('Insert size')
  plt.ylabel('Frequency')
  plt.title( 'Histogram of insert size for %s' % fname )
  plt.text(mu+20, .025, 'mean=%.2f\nstdev=%.2f\nmedian=%s' % (mu,sd,median) ) # plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
  plt.axis([0, 2*mu, 0, 0.04])
  #plt.grid(True)
  #plt.show()
  plt.savefig( fname+'.png',dpi=300 )

def main( args ):
  fname=args[0]
  maxInsert=999999
  if len(args)>1:
    maxInsert=int(args[1])
  if len(args)>1: maxInsert=int(args[1])
  
  handle=sys.stdin
  insertsSizes=[]
  for line in handle:
    QNAME,FLAG,RNAME,POS,MAPQ,CIAGR,MRNM,MPOS,TLEN,SEQ,QUAL=line.split('\t')[:11]
    TLEN=fabs( int(TLEN) )#; print TLEN
    if TLEN and TLEN<maxInsert: insertsSizes.append( TLEN )
    
  histplot( insertsSizes,fname )
  
if __name__=='__main__': 
  t0=datetime.now()
  main( sys.argv[1:] )
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
