#!/usr/bin/env python2
"""

Author:
l.p.pryszcz+git@gmail.com

Barcelona, 22/04/2012
"""
 
import os, sys
from datetime import datetime
from optparse import OptionParser

def main():
  
    usage = "usage: cat bigfile | %prog [options]" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True
  
    parser.add_option("-f", dest="files", default=2, type=int,
                      help="split into f files  [%default]" )
    parser.add_option("-l", dest="lines", default=4, type=int,
                      help="split every l lines [%default]" )
    parser.add_option("-n", dest="name",  default="out", type=str,
                      help="output name         [%default]" )
    parser.add_option("-e", dest="ext",   default=".fq", type=str,
                      help="output extension    [%default]" )
  
    ( o, fPaths ) = parser.parse_args()
    files,lines = o.files,o.lines
    
    #open out files    
    outFiles=[]
    for i in range( 1,files+1 ):
        outFiles.append( open("%s%s%s" % (o.name,i,o.ext),"w") )

    #split lines
    li=fi=0
    for l in sys.stdin:
       if li==lines:
           li=0
           fi+=1
           if fi==files:
               fi=0
       outFiles[fi].write( l )
       li += 1

    #close files
    for f in outFiles:
        f.close()

if __name__=='__main__': 
    t0=datetime.now()
    main()
    dt=datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
