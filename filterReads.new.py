#!/usr/bin/env python
desc="""Filter QSEQ/FastQ reads. Store output as FastQ.
Reads are clipped at first undetermined base (. in qseq or N in fastq)
and at first base having qual below -q.
Reads (and their pairs if -p) not passing filtering are discarded. 
Orphaned reads may be store optionally (-u). 
"""
epilog="""Author:
Leszek Pryszcz
l.p.pryszcz@gmail.com

Barcelona, 5/11/2013
"""

"""
Fixes:
-0.2
--preformance:
---multiprocessing tested, but very little improvement
---zcat subprocess speeds up by 30% (1M of paired 100bp reads processed in 40s in neptune)
--tested biopython, but ~4x slower than raw fq parser

-0.1
--wrong sep in _clipSeq solved

-0.11
--output always PHRED+33 quals (Sanger, CASAVA1.8+)
--include reads with '.' bases -> 'N'
"""

import argparse, gzip, os, sys
import subprocess
from datetime import datetime
from Bio import SeqIO
import locale
locale.setlocale(locale.LC_ALL, 'en_US.utf8')

def qseqparser(handle):
    """Parse QSEQ fromat and yield name, sequence and qualities."""
    for l in handle:
        qseq_element = l[:-1].split('\t') #SOLEXA 90403 4 1 23 1566 0 1 ACCGCTCTCGTGCTCGTCGCTGCGTTGAGGCTTGCG `aaaaa```aZa^`]a``a``a]a^`a\Y^`^^]V` 1
        if len(qseq_element) != 11 or qseq_element[-1] != '1': 
            yield
            continue
        #formatting
        name       = '@%s:%s:%s:%s:%s#%s/%s' % (qseq_element[0], qseq_element[2], qseq_element[3], qseq_element[4], qseq_element[5], qseq_element[6], qseq_element[7])
        seq, quals = qseq_element[8], qseq_element[9]
        seq        = seq.replace(".", "N")
        yield name, seq, quals

def fqparser(handle):
    """Parse fastq format and yield name, sequence and qualities."""
    fqlist = []
    for l in handle:
        fqlist.append(l[:-1])
        if len(fqlist) != 4:
            continue
        #unload & reset fqlist
        name, seq, sep, quals = fqlist
        fqlist = []
        yield name, seq, quals
        
def _clipSeq(seq, quals, sep='.'):
    """Clip sequence at first sep base (. or N). Clip quals accordingly."""
    if sep in seq:
        pos = seq.index(sep)
        seq, quals = seq[:pos], quals[:pos]
    return seq, quals
        
def rawtrimmer(infile, minlen, minqual, qual64offset, qseq, stripHeaders, outformat, pi, pair=""):
    """Single process implementation of rawtrimmer.
    Open zcat subprocess and read from stdin."""
    handle = infile
    if infile.name.endswith('.gz'):
        zcat   = subprocess.Popen(['zcat', infile.name], bufsize=-1, stdout=subprocess.PIPE)
        handle = zcat.stdout
        #handle = gzip.open(infile.name)
    #get parser
    if qseq:
        parser = qseqparser(handle)
    else:
        parser = fqparser(handle)
    #process entries
    for seqi, (name, seq, quals) in enumerate(parser, pi+1):
        #clip seq & quals @ N ( unknown base )
        seq, quals = _clipSeq(seq, quals, 'N')
        #check if correct length
        if not seq or len(seq) < minlen:
            yield
            continue
        #return PHRED+33 quals (Sanger encoding)
        if qual64offset:
            quals=''.join([chr(ord(q)-31) for q in quals])
        #cut sequence & quals @ quality
        if minqual:
            for pos, qual in enumerate(quals):
                if ord(qual)-33 < minqual:
                    seq = seq[:pos]
                    break
            #check if correct length
            if not seq or len(seq) < minlen:
                yield
                continue
            #clip seq and qual
            quals = quals[:len(seq)]
        #define fastQ line
        if stripHeaders:
            name = "@%s%s" % (seqi, pair)
        if outformat == "fasta":
            fastq = '>%s\n%s\n' % (name[1:], seq)
        else:
            fastq = '%s\n%s\n+\n%s\n' % (name, seq, quals)
        yield fastq
    
def filter_paired(fpair, outfiles, minlen, minqual, qual64offset, qseq, \
                  stripHeaders, outformat, pi):
    """Filter paired reads."""
    inF, inR = fpair
    outF, outR, outCombined, outUnpaired = outfiles
    #define parsers rawtrimmer fqtrimmer
    fqparser1 = rawtrimmer(inF, minlen, minqual, qual64offset, qseq, stripHeaders, outformat, pi)
    fqparser2 = rawtrimmer(inR, minlen, minqual, qual64offset, qseq, stripHeaders, outformat, pi)
    #process 
    both = fori = revi = filtered = 0
    for i, rec1 in enumerate(fqparser1, pi+1):
        #will crash if len(fq1) > len(fq2)
        rec2 = fqparser2.next()
        if rec1 and rec2:
            #store paired output
            if outF:
                outF.write(rec1)
                outR.write(rec2)
            #store combined output
            if outCombined: 
                outCombined.write(rec1+rec2)
            both += 1
        elif outUnpaired and rec1:
            #store F read if R didn't pass filtering and orphans requested
            fori+=1
            outUnpaired.write(rec1)
        elif outUnpaired and rec2:
            #store R read if F didn't pass filtering and orphans requested
            revi+=1
            outUnpaired.write(rec2)
        else:
            #count skipped reads
            filtered += 1             
        #print stats
        if not i % 10e3:
            info = "%9s processed [%6.2f%s ok] Both passed: %s Orphans F/R: %s/%s      \r" % (locale.format("%d", i, grouping=True), (i-filtered)*100.0/i, '%', both, fori, revi)
            sys.stderr.write(info)
    info = "%9s processed [%6.2f%s ok] Both passed: %s Orphans F/R: %s/%s      \n" % (locale.format("%d", i, grouping=True), (i-filtered)*100.0/i, '%', both, fori, revi)
    sys.stderr.write(info)
    return i, filtered, fori+revi

def process_paired(inputs, qseq, outdir, unpaired, minlen, minqual, \
                   noSeparate, combined, qual64offset, replace, \
                   stripHeaders, fasta, verbose):
    """Process paired libraries."""
    #define output fnames
    fnend = outformat = 'fastq'
    if fasta:
        fnend = outformat = 'fasta'
    outfnF     = os.path.join(outdir, 'q%s_1.%s'        % (minqual, fnend))
    outfnR     = os.path.join(outdir, 'q%s_2.%s'        % (minqual, fnend))
    unpairedfn = os.path.join(outdir, 'q%s.unpaired.%s' % (minqual, fnend))
    combinedfn = os.path.join(outdir, 'q%s.combined.%s' % (minqual, fnend))      
    #check if outfiles exists
    if not replace:
        if os.path.isfile(outfnF) or os.path.isfile(outfnR) or os.path.isfile(unpairedfn) or os.path.isfile(combinedfn):
            sys.stderr.write("At least one of the output files is present. Remove them or run with --replace parameter. Exiting!\n")
            sys.exit()
            
    #open files for writting
    outF = outR = outCombined = outUnpaired = False
    if not noSeparate:
        outF = open(outfnF, 'w')
        outR = open(outfnR, 'w')
    #open out file for unpaired reads
    if unpaired: 
        outUnpaired = open(unpairedfn, 'w')
    #open out file for combined FastQ
    if combined: 
        outCombined = open(combinedfn, 'w')
    outfiles = (outF, outR, outCombined, outUnpaired)
    
    #process all input files
    fpair = []
    i = pi = filtered = single = 0
    for fn in inputs:
        fpair.append(fn)
        if len(fpair) != 2:
            continue
        #proces qseq files: GERALD->FASTA
        i, pfiltered, psingle = filter_paired(fpair, outfiles, minlen, minqual, qual64offset, qseq, stripHeaders, outformat, pi)
        #print info
        if verbose: 
            sys.stdout.write('[%s]  %s  %s  %s  %s\n' % (datetime.ctime(datetime.now()), fpair[0].name, fpair[1].name, i-pi, pfiltered))
        #update read counts
        pi        = i
        filtered += pfiltered
        single   += psingle
        #reset fnames
        fpair = []
      
    #close outfiles
    for outfile in outfiles:
        if outfile:
            outfile.close()
    #print info 
    sys.stdout.write('Processed pairs: %s. Filtered: %s. Reads pairs included: %s [%.2f%s]. Orphans: %s [%.2f%s]\n' % ( i,filtered, i-filtered, (i-filtered)*100.0/i, '%', single, single*100.0/i,'%'))

def filter_single(infile, out, minlen, minqual, qual64offset, qseq, \
                  stripHeaders, outformat, pi):
    """Filter single reads."""
    #define parser
    fqparser = fqtrimmer(infile, minlen, minqual, qual64offset, qseq, stripHeaders, outformat, pi)
    #process output of subprocesses
    both = filtered = 0
    for i, rec in enumerate(fqparser, pi+1):
        #write pair if F & R passed filtering
        if rec:
            out.write(rec)
            both += 1
        #nothing if both didn't pass filtering 
        else: 
            filtered += 1             
        #print stats
        if not i % 10e3:
            info = "%9s processed [%6.2f%s ok]      \r" % (locale.format("%d", i, grouping=True), (i-filtered)*100.0/i, '%')
            sys.stderr.write(info)
    info = "%9s processed [%6.2f%s ok]       \n" % (locale.format("%d", i, grouping=True), (i-filtered)*100.0/i, '%')
    sys.stderr.write(info)
    return i, filtered
    
def process_single(inputs, qseq, outdir, minlen, minqual, qual64offset, \
                   replace, stripHeaders, fasta, verbose):
    """Process single end libraries."""
    #define output fnames
    fnend = outformat = 'fastq'
    if fasta:
        fnend = outformat = 'fasta'
    outfn     = os.path.join(outdir, 'q%s.%s' % (minqual, fnend))
    #check if outfiles exists
    if not replace:
        if os.path.isfile(outfn):
            sys.stderr.write("File exists: %s. Remove them or run with --replace parameter. Exiting!\n"%outfn)
            sys.exit()
    #process input files
    i = pi = filtered = 0
    out = open(outfn, 'w')
    for fn in inputs:
        #proces qseq files: GERALD->FASTA
        i, pfiltered = filter_single(fn, out, minlen, minqual, qual64offset, qseq, stripHeaders, outformat, pi)
        #print info
        if verbose: 
            sys.stdout.write('[%s]   %s  %s  %s\n' % (datetime.ctime(datetime.now()), fn, i-pi, pfiltered))
        #update read counts
        pi        = i
        filtered += pfiltered
    #close outfile
    out.close()
    #print info 
    sys.stdout.write('Processed: %s. Filtered: %s. Reads included: %s [%.2f%s].\n' % (i, filtered, i-filtered, (i-filtered)*100.0/i, '%'))        
    
def main():
  
    usage  = "%(prog)s -p -u -q10 -l31 -i sample_read1.fastq.gz sample_read2.fastq.gz -o sample [options]" 
    parser = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", "--verbose", default=False, action="store_true")
    parser.add_argument("--version", action="version", default="0.2")
    parser.add_argument("-i", "--inputs", nargs="+", type=file,
                        help="input file(s)")
    parser.add_argument("-o", "--outdir", default='outdir', 
                        help="define where to store output files")
    parser.add_argument("-g", "--qseq", action="store_true",  default=False,
                        help="input is QSEQ, not FastQ")
    parser.add_argument("-l", "--minlen", default=31, type=int,
                        help="min read lenght (shorter after quality trimming are removed) [%(default)s]" )
    parser.add_argument("-q", "--minqual", default=None, type=int,
                        help="read is clipped @ first base having PHRED quality lower than [%(default)s]" )
    parser.add_argument("-t", dest="qual64offset", default=False, action='store_true',
                        help="use PHRED+64 (illumina/solexa) quality encoding [Sanger PHRED+33 encoding]")
    parser.add_argument("-p", "--paired", default=False, action="store_true", 
                        help="paired-end reads (qXX_1.fastq & qXX_2.fastq)")
    parser.add_argument("-u", "--unpaired", default=False, action="store_true", 
                        help="store orphaned reads > qXX.unpaired.fastq")                  
    parser.add_argument("-r", "--replace", default=False, action="store_true", 
                        help="overwrite output files")
    parser.add_argument("-b", "--noSeparate", default=False, action="store_true", 
                        help="don't store separate fastQ for F & R reads > qXX_1.fastq qXX_2.fastq" )                  
    parser.add_argument("-c", "--combined", default=False, action="store_true", 
                        help="store combined fastQ for paired reads > qXX.combined.fastq [%(default)s]" )
    parser.add_argument("-H", "--stripHeaders", default=False, action="store_true", 
                        help="replace headers by int [%(default)s]" )
    parser.add_argument("--fasta", default=False, action="store_true", 
                        help="report fasta, not FastQ" )
  
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n" % str(o))
    
    #create output directore if not present already
    if not os.path.isdir(o.outdir): 
        os.makedirs(o.outdir)

    #check if quality encoding correct
    ##Quality offset checking need to be done!
        
    #gerald2fastq
    if o.paired:
        process_paired(o.inputs, o.qseq, o.outdir, o.unpaired, o.minlen, o.minqual, \
                       o.noSeparate, o.combined, o.qual64offset, o.replace, \
                       o.stripHeaders, o.fasta, o.verbose)
    else:
        process_single(o.inputs, o.qseq, o.outdir, o.minlen, o.minqual, \
                       o.qual64offset, o.replace, \
                       o.stripHeaders, o.fasta, o.verbose )
                  
if __name__=='__main__': 
    t0=datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!                \n")
    dt=datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
