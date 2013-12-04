#!/usr/bin/env python
"""Parses multiple files and find common rows between these.

Author:
l.p.pryszcz@gmail.com

Barcelona, 21/03/2012
"""

import os, sys
from optparse import OptionParser
from datetime import datetime

def compare_csv( fnames,exclude,include,fraction,columns,delimiter,minCommon,splitFn ):
    """
    """
    exclude = set( exclude )
    include = set( include )
    #define columns
    if columns.isdigit(): #select column range if digit
        columns = int( columns )
    else:                 #otherwise select given number of columns
        columns = [ int(x) for x in columns.split(',') ]
    #process files
    common = {}
    for fn in fnames:
        fnShort = fn
        if splitFn:
            fnShort = fn.split(".")[0]

        for l in open(fn):
            l = l.strip()
            if l.startswith("#") or not l:
                continue
            
            lData = l.split(delimiter)
            if type(columns)==int:
                k = tuple( lData[:columns] )
            else:
                k = []
                for c in columns:
                    if c >= len(lData):
                        continue
                    k.append( lData[c] )
                k = tuple( k )

            if k not in common:
                common[k] = set()
            
            common[k].add( fnShort )
    
    c2i = {}
    print "#\tsamples\tfields"
    for k in sorted( common.keys(),key=lambda x: len(common[x]), reverse=True ):
        c = len(common[k])
        if c not in c2i:
            c2i[c] = 0
        c2i[c] += 1
        # break if less common than defined
        if c<minCommon:
            continue # break
        
        # check exclude/include info
        if exclude and include:
            # get fraction of include and exclude
            inFract = len(include.intersection(set(common[k])))*1.0/len(include)
            exFract = len(exclude.intersection(set(common[k])))*1.0/len(exclude)
            
            if   inFract>=fraction   and exFract<=1-fraction: 
                pass # enough include and exclude less than 1-fraction
            elif inFract<=1-fraction and exFract>=fraction: 
                pass # include less than 1-fraction and enough exclude
            else:
                continue
                
        # print info
        print "%s\t%s\t%s" % ( c,",".join( common[k] ),delimiter.join( k ) )

    sys.stderr.write( "#common for\t#\n" )
    for c in sorted( c2i.keys(),reverse=True ):
        sys.stderr.write( "%s\t%s\n" % (c,c2i[c] ) )

def main():
    
    usage = "usage: %prog [options]" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True

    parser.add_option("-e", dest="exclude", default="", 
                      help="exclude")
    parser.add_option("-i", dest="include", default="", 
                      help="include")
    parser.add_option("-f", dest="fraction", default=1.0, type=float,
                      help="exclude/include fraction [%default]")
    parser.add_option("-c", dest="columns",   default=1, type=str, 
                      help="""if int: this number of first column will be compared;
 or comma separated column numbers (0-based) [%default]""")
    parser.add_option("-d", dest="delimiter", default="\t", 
                      help="column delimiter [tab]")
    parser.add_option("-m", dest="minCommon", default=1, type=int,
                      help="# repeats to be reported [%default]")
    parser.add_option("-s", dest="splitFn",  default=False, action="store_true", 
                      help="split fname (sheet name) by dot")
    parser.add_option("-v", dest="verbose",  default=True, action="store_false")
    
    ( o, args ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "%s\n" % ( str(o), ) )

    # check if files exists
    for fn in args:
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" )

    # define exclude/include
    if o.exclude and o.include:
        exclude = o.exclude.split(',')
        include = o.include.split(',')
    else:
        exclude,include = [],[]

    # process csv files
    compare_csv( args,exclude,include,o.fraction,o.columns,o.delimiter,o.minCommon,o.splitFn )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
