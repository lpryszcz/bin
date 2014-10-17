#!/usr/bin/env python
desc="""Report random entries from FastQ/FastA file(s).
Support uncompressed or BGZIP compressed files. 
Headers need to have the same name in FastQ/FastA files. 
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 13/10/2014
"""

import gzip, os, random, sqlite3, sys, time
from datetime import datetime
from Bio import SeqIO, bgzf

class Counter(object):
    """Counter class"""
    def __init__(self, count=-1):
        self.reset_counter(count)
        
    def reset_counter(self, count=0):
        """Reset counter to given count."""
        self.counter = count

    def count(self, *args):
        """Return function calls"""
        self.counter += 1
        return self.counter

def get_seq(handle, offset, length):
    """Return entry for given file, offset and length"""
    handle.seek(offset)
    return handle.read(length)

def store_random_entries(outbase, infiles, n, verbose):
    """Return target fastas for protids from sqlite3."""
    if verbose:
        sys.stderr.write("Preparing files...\n")
    files = []
    cursors = []
    outfiles = []
    #open target files/cursors
    for fi, f in enumerate(infiles, 1):
        #get & store cursor
        cur = sqlite3.connect(f.name+'.idx').cursor()
        cursors.append(cur)
        #get files
        cur.execute("SELECT name FROM file_data")
        files.append({})
        #open outfiles
        outfiles.append(gzip.open(outbase+".%s.fq.gz"%fi, "w"))
        for name, in cur.fetchall():
            #fpath = os.path.join(os.path.dirname(db), name)
            if name.endswith('.gz'):
                files[-1][name] = bgzf.open(name)
            else:
                files[-1][name] =      open(name)

    #preload offset data for other files
    if verbose:
        sys.stderr.write("Loading offset_data...\n")                
    cmd1 = """SELECT f.name, offset, length FROM offset_data o JOIN file_data f
    ON o.file_number=f.file_number"""
    offset_data = [[], ]
    for cur in cursors[1:]:
        cur.execute(cmd1)
        offset_data.append(cur.fetchall())
    #get randomly sorted first file
    if verbose:
        sys.stderr.write("Selecting random entries...\n")
    cmd0 = """SELECT key, f.name, offset, length FROM offset_data o JOIN file_data f
    ON o.file_number=f.file_number ORDER BY RANDOM()"""
    if n>0:
        cmd0 += " LIMIT %s"%n
    cursors[0].execute(cmd0)
    #combine randomised and preloaded data
    if verbose:
        sys.stderr.write("Reporting...\n")    
    for i, (key, name, offset, length) in enumerate(cursors[0].fetchall(), 1):
        if verbose and i%10000==1:
            sys.stderr.write(" %s     \r"%i)
        #store first file random sequence
        outfiles[0].write(get_seq(files[0][name], offset, length))
        #store sequence from the remaining files
        key = int(key)
        for fi in range(1, len(outfiles)):
            name, offset, length = offset_data[fi][key]
            outfiles[fi].write(get_seq(files[fi][name], offset, length))
    #close
    for out in outfiles:
        out.close()
        
def fastq2random(outbase, files, n, verbose, seqformat='fastq'):
    """Return number of random reads from FastQ file(s)"""
    #generate indexes
    if verbose:
        sys.stderr.write("Generating indexes...\n")
    indexes = [] 
    for i, f in enumerate(files, 1):
        sys.stderr.write(" %s \r"%i)
        c = Counter() #counters.append(Counter())
        index = SeqIO.index_db(f.name+'.idx', f.name, seqformat, key_function=c.count)
        index.close()

    #get random entries
    store_random_entries(outbase, files, n, verbose)

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "--input", nargs="+", type=file, 
                        help="input file(s)")
    parser.add_argument("-o", "--outbase", 
                        help="output base file name")
    parser.add_argument("-n", type=int, default=0, 
                        help="report n entries [all in randomised order]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    fastq2random(o.outbase, o.input, o.n, o.verbose)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
