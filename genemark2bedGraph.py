#!/usr/bin/env python2
"""Convert gdata from gene mark into bedGraph.

USAGE:
genemark2bedGraph.py *.gdata [> bedGraph]

"""

import os,sys

def fn2bedGraph( fn ):
  """Generate bedGraph
  http://genome.ucsc.edu/goldenPath/help/bedgraph.html
  """
  #chromosome name is chr1.gdata -> chr1
  dotIndex = fn.rindex(".")
  s = fn[:dotIndex]
  contig=s
  #iterate through windows
  start=0
  for l in open(fn):
    #skip empty and commented lines
    l = l.strip()
    if not l or l.startswith("#"):
      continue
    #get data from line
    lData = [ float(x) for x in l.split("\t") ]
    end   = int( lData[0] )
    #take max probability
    prob  = max( lData[1:] )
    #print output
    print "%s\t%s\t%s\t%.3f" % ( contig,start,end,prob )
    #define new start
    start = end

if __name__=="__main__":
  fnames = sys.argv[1:]
  if not fnames:
    sys.exit( "USAGE:\ngenemark2bedGraph.py chr1.gdata [ chr2.gdata > bedGraph]" )
  
  #iterate through files
  print "track type=bedGraph name='coding potential' graphType=bar viewLimits=0:1 windowingFunction=mean"
  for fn in fnames:
    if not os.path.isfile( fn ):
      sys.exit( "No such file: %s" % fn )
    sys.stderr.write( " %s       \r" % fn )
    #generate bedGraph
    fn2bedGraph( fn )
    
