#!/usr/bin/env python
"""Generate Latex template based on xls file. 

Author:
l.p.pryszcz@gmail.com

Dolny Smokovec, 19/08/2016
"""

import os, sys
from optparse import OptionParser
from datetime import datetime

try:
    import xlrd
except:
    sys.stderr.write("Install xlrd first: sudo apt-get install python-xlrd\n")
    sys.exit(1)

try:
    from latexencode import utf8tolatex
except:
    sys.stderr.write("Install pylatexenc first: sudo pip install pylatexenc && cp /usr/local/lib/python2.7/dist-packages/pylatexenc/latexencode.py .\n")
    sys.exit(1)

def xls2tex(xls, out):    
    """Convert xls to tex"""
    # define empty workbook
    book = xlrd.open_workbook(xls)
    sh   = book.sheet_by_index(0)

    tex = u"""\\title{%s}{%s}\n{%s}{%s}{%s}{%s}\n\n"""
    section = "Abstracts"
    time = ""
    for rx in range(sh.nrows):
        if sh.row(rx)[0].value.startswith('#'):
            continue
        cells = [utf8tolatex(c.value) for c in sh.row(rx)]
        author, email, affiliation, title, abstract = cells
        out.write(tex%(section, title, author, affiliation, time, abstract))
    
def main():
    
    usage = "usage: %prog [options]" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True

    parser.add_option("-i", "--xls",     default="abstracts.xls",
                      help="output fname     [%default]")
    parser.add_option("-o", "--out",     default="abstracts.tex",
                      help="output fname     [%default]")
    parser.add_option("-v", dest="verbose",  default=True, action="store_false")
    
    ( o, args ) = parser.parse_args()

    with open(o.out, "w") as out:
        xls2tex(o.xls, out)
        
if __name__=='__main__': 
  t0=datetime.now()
  main() 
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
