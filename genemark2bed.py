#!/usr/bin/env python2
"""Convert gdata from gene mark into bedGraph.

USAGE:
genemark2bedGraph.py *.gdata [> bedGraph]

AUTHOR:
l.p.pryszcz+git@gmail.com

Barcelona, 17/07/2012
"""

import os,sys

def fn2bed( fn ):
  """Generate bed
  http://genome.ucsc.edu/goldenPath/help/bedgraph.html
  """
  #chromosome name is chr1.gdata -> chr1
  dotIndex = fn.rindex(".")
  s = fn[:dotIndex]
  contig=s
  #iterate through windows
  for l in open(fn):
    #skip empty and commented lines
    l = l.strip()
    if not l or l.startswith("#"):
      continue
    lData = l.split()
    if len(lData)!=5:
      continue
    #get data from line
    start,end,frame,prob,startProb = lData
    try:
      start,end,frame = int(start),int(end),int(frame)
      prob,startProb  = float(prob),float(startProb)
    except:
      sys.stderr.write( "Wrong line: %s\n" % lData )
      continue
    #get frame
    #  Frames are counted 1-6, Direct 1-3 and Complement 4-6
    #  The first base of the input sequence is numbered 1 not 0
    #  The frame is computed by the direct strand gene start 
    #    (LEnd for direct strand, REnd for complement strand)
    #     Frame=(LEnd-1)%3+1  or (REnd-1)%3+4
    strand = "+"
    if frame>3:
      strand = "-"
    #print output
    print "%s\t%s\t%s\tsp=%s\t%.3f\t%s" % ( contig,start-1,end,startProb,prob,strand )

if __name__=="__main__":
  fnames = sys.argv[1:]
  if not fnames:
    sys.exit( "USAGE:\ngenemark2bed.py chr1.ldata [ chr2.ldata > genemark.bed ]" )
  
  #iterate through files
  print "track type=bed name=genemark"
  for fn in fnames:
    if not os.path.isfile( fn ):
      sys.exit( "No such file: %s" % fn )
    sys.stderr.write( " %s       \r" % fn )
    #generate bedGraph
    fn2bed( fn )
    
