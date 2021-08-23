#!/usr/bin/env python2
# Convert arabic chr (IV) into int (4).
# Can handle chr. 

import sys
from functools import reduce

numerals = {'M' : 1000, 'D' : 500, 'C' : 100, 'L' : 50, 'X' : 10, 'V' : 5, 'I' : 1}
def romannumeral2number(s):
    return reduce(lambda x, y: -x + y if x < y else x + y, map(lambda x: numerals.get(x, 0), s.upper()))

handle = sys.stdin
out = sys.stdout
for l in handle:
    #if l.startswith(">"):
    ref = l.split()[0]
    addchr = ""
    if ref.lower().startswith('chr'):
        addchr = ref[:3]
        ref = ref[3:]
    nref = romannumeral2number(ref)
    if addchr: ref = addchr+ref
    l = "%s"%nref + l[len(ref):]
    out.write(l)
