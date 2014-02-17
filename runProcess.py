#!/usr/bin/env python
###
# Get parameters and run given command.
#
###
import os
import sys

def runProcess( argv ):
  #print argv
  command=''
  for arg in argv:
    command+=arg+' '
  #print command
  os.system(command)
  
if __name__=='__main__': sys.exit(runProcess(sys.argv[1:]))
