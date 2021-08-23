#!/usr/bin/env python2
"""
xlsx2tsv file1.xlsx [file2.xlsx ... fileN.xlsx]

Parse a .xlsx file(s) (Excel OOXML, which is not OpenOffice) into tab-separated values.
Output directory (ie. file1.xls_sheets/) is created per each input file.
In this directory the sheets from given xlsx file are stored as separate .tsv files. 
Outputs honest-to-goodness tsv, no quoting or embedded \\n\\r\\t.

The spec for this format is 5220 pages.  I did not use it.  This was helpful:
http://blogs.msdn.com/excel/archive/2008/08/14/reading-excel-files-from-linux.aspx
But mostly I guessed how the format works.  So bugs are guaranteed.

# by brendan o'connor - anyall.org - gist.github.com/22764
# Modified by lpryszcz (http://bioinfoexpert.com/) 
"""

import xml.etree.ElementTree as ET
import os, sys, zipfile, re, itertools

if len(sys.argv)<2:
  print __doc__.strip()
  sys.exit(1)

def letter2col_index(letter):
  """ A -> 0, B -> 1, Z -> 25, AA -> 26, BA -> 52 """
  base26digits = [1+ord(x)-ord("A") for x in letter]
  return sum([x*26**(len(base26digits) - k - 1)  for k,x in enumerate(base26digits)]) - 1

def flatten(iter):
  return list(itertools.chain(*iter))

def cell2text(cell, ss_list):
  if cell is None:
    return ""
  elif 't' in cell.attrib and cell.attrib['t'] == 's':
    # shared string
    idx = int(cell.find(n("v")).text)
    si = ss_list[idx]
    t_elt = si.find(n("t"))
    if t_elt is not None:
      return t_elt.text
    t_elts = si.findall(n("r") + "/" + n("t"))
    if t_elts:
      text = "".join( (t.text) for t in t_elts )
      return text
    raise Exception("COULDNT DECODE CELL: %s" % ET.tostring(si))
    #return si.find(n("t")).text
    #return ET.tostring(si)
  else:
    v_elt = cell.find(n("v"))
    if v_elt is None: return ""
    return v_elt.text

def make_cells(max_col):
  return [None] * (max_col+1)

warning_count=0
warning_max = 50
def warning(s):
  global warning_count
  warning_count += 1
  if warning_count > warning_max: return
  print>>sys.stderr, "WARNING: %s" % s

def cell_text_clean(text):
  s = text.encode("utf-8")
  if "\t" in s: warning("Clobbering embedded tab")
  if "\n" in s: warning("Clobbering embedded newline")
  if "\r" in s: warning("Clobbering embedded carriage return")
  s = s.replace("\t"," ").replace("\n"," ").replace("\r"," ")
  return s

def read_sheet_and_save(z, sheet_num, out):
  ss_xml = z.read("xl/sharedStrings.xml")
  ss_list = ET.XML(ss_xml).findall(n("si"))

  xml = z.read("xl/worksheets/sheet%s.xml" % sheet_num)
  s = ET.fromstring(xml)
  rows = s.findall(n("sheetData")+"/"+n("row"))

  all_cells = flatten( [[c for c in row.findall(n("c"))] for row in rows] )
  max_col = max(letter2col_index(re.search("^[A-Z]+",c.attrib['r']).group()) for c in all_cells)

  
  for row in rows:
    cells_elts = row.findall(n("c"))
    inds = []  # parallel
    for c in cells_elts:
      letter = re.search("^[A-Z]+", c.attrib['r']).group()
      inds.append(letter2col_index(letter) )
    cells = make_cells(max_col)
    for c,j in zip(cells_elts,inds): 
      cells[j] = c
    #print( *(cell2text( c ).encode("utf-8").replace("\t"," ") for c in cells), sep="\t")
    #print myjoin((cell_text_clean(cell2text( c )) for c in cells), sep="\t")
    out.write("\t".join(cell_text_clean(cell2text(c, ss_list)) for c in cells)+"\n")  

  if warning_count > warning_max:
    print>>sys.stderr, "%d total warnings, %d hidden" % (warning_count, warning_count-warning_max)

# define lambda
n=lambda x: "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}%s" % x
    
# parse multiple xlsx files
sys.stderr.write("Processing xlsx file(s)...\n")
for fname in sys.argv[1:]:
  sys.stderr.write(" %s\n"%fname)
  # prepare output directory
  outdir = fname+"_sheets"
  if not os.path.isdir(outdir):
    os.makedirs(outdir)
  # load file
  z = zipfile.ZipFile(fname)
  # get sheets
  sheet_filenames = [f for f in z.namelist() if re.search("^xl/worksheets/sheet.*xml$", f)]
  workbook_x = ET.XML(z.read("xl/workbook.xml"))
  sheet_xs = workbook_x.find(n("sheets")).findall(n("sheet"))
  # save sheets
  for sheet_num, x in enumerate(sheet_xs, 1):
    name = x.get('name')
    outfn = os.path.join(outdir, name)+'.tsv'
    sys.stderr.write("  %s %s > %s\n"%(sheet_num, name, outfn))
    out = open(outfn, "w") #sys.stdout
    read_sheet_and_save(z, sheet_num, out)
