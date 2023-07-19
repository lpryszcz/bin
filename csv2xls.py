#!/usr/bin/env python3
desc="""Convert multiple cvs files into xls.
Each csv file becomes one sheet in xls.
By default, tab-delimited (tabbed) files are processed. But it 
can be change with -d parameter.

Pre-requisities:
xlwt (http://scienceoss.com/write-excel-files-with-python-using-xlwt/)

Note: As Excell 97-2003 stores row number as 16 bit integer, 
it can handle 65536 rows at max! There is not way around this:/

"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 21/03/2012
"""

import os, sys
#from optparse import OptionParser
from datetime import datetime

#try to import module for xls handling
try:
    from xlwt import * #import xlwt #http://scienceoss.com/write-excel-files-with-python-using-xlwt/
except:
    #sys.stderr.write("xlwt not installed!\nHave a look at how to install it on: http://scienceoss.com/write-excel-files-with-python-using-xlwt/\n")
    sys.stderr.write("Install xlwt first: sudo apt-get install python-xlwt OR conda install xlwt\n")
    sys.exit(1)
##Excell formatting
f = Font()
#f.height = 20*72
#f.name = 'Verdana'
#f.bold = True
f.italic = True
#f.underline = Font.UNDERLINE_DOUBLE
f.colour_index = 4

h_style = XFStyle()
h_style.font = f


def format_value( s ):
    """Convert str to int or float.
    """
    #digit?
    if s.isdigit():
        return( int(s) )
    #float?
    try:
        return float(s)
    except:
        pass
    #string...
    return s

def add_sheet( wbk,fn,delimiter,splitFn ):
    """
    """
    sname = fn.replace("/","_").replace("-","")
    if splitFn:
        sname = sname.split('.')[0]
    #add new sheet
    sheet = wbk.add_sheet( sname,cell_overwrite_ok=False )
    #add lines
    row = 0
    for l in open( fn ):
        col = 0
        #add columns to row
        for val in l[:-1].split( delimiter ):
            val = format_value( val )
            sheet.write(row,col,val)            
            col += 1
        row += 1

    return wbk

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-i", "--fnames", nargs="+", 
                      help="input csv files")
    parser.add_argument("-o", "--outfn",     default="out.xls",
                      help="output fname     [%(default)s]")
    parser.add_argument("-d", "--delimiter", default="\t", 
                      help="column delimiter [tab]")
    parser.add_argument("-s", "--splitFn",  action="store_true", 
                      help="split fname (sheet name) by dot")
    parser.add_argument("-v", "--verbose", action="store_false")
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))
    #quit if out file exists
    if os.path.isfile( o.outfn ):
        parser.error( "File exists: %s" % o.outfn )
    #check if files exists
    for fn in o.fnames:
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" )
        #check if no more than 65536 lines - xls cannot handle more!
        i=0
        for l in open(fn):
            i+=1
            if i>65536:
                parser.error( "More than 65536 in %s" % fn )

    #define empty workbook
    wbk = Workbook()
    #add all files to workbook
    for fn in o.fnames:
        wbk = add_sheet( wbk,fn,o.delimiter,o.splitFn )
    #write workbook to disc
    wbk.save( o.outfn )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
