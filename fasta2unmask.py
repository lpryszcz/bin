#!/usr/bin/env python3
import os, sys

out = sys.stdout
for l in sys.stdin:
    if l.startswith(">"): out.write(l)
    else: out.write(l.upper())
    
