#!/usr/bin/env python
"""
Identify SNP sites in mpileup out from BAM alignments. 
In addition calculates overall and each contig coverage. 
 
At first run: /users/tg/lpryszcz/cluster/assembly/mpileup.sh (s=MCO456; samtools mpileup -f gem/$s.fa gem/${s}.bam > gem/${s}.mpileup)
Execute: 
 python ~/workspace/assembly/src/heterozygosity.py -i gem/CBS6318.mpileup [-f minFreq -d minDepth -o out_base_name]
or use with piping:
 s=CBS1954; samtools mpileup -f gem/$s.fa gem/$s.bam | python ~/workspace/assembly/src/heterozygosity.py -o gem/$s -f 0.2 -d 10
"""

import os, sys
from optparse import OptionParser
from datetime import datetime
import subprocess

def _remove_indels( alts ):
    """
    Remove indels from mpileup.
    .$....,,,,....,.,,..,,.,.,,,,,,,....,.,...,.,.,....,,,........,.A.,...,,......^0.^+.^$.^0.^8.^F.^].^],
    ........,.-25ATCTGGTGGTTGGGATGTTGCCGCT..
    """
    #remove indels info
    for symbol in ('-','+'):
        baseNo = 0
        while symbol in alts:
            i=alts.index(symbol)
      
            j = 1
            digits=[]
            while alts[i+j].isdigit():
                digits.append( alts[i+j] )
                j += 1
      
            if digits:
                baseNo=int( ''.join(digits) )
        
            alts=alts[:i]+alts[i+baseNo+len(digits)+1:] #......+1A..,
      
    return alts

def get_alt_allele( base_ref,cov,alg,minFreq,alphabet,reference,bothStrands ):
    """Return alternative allele only if different than ref and freq >= minFreq.
    """
    #remove deletions
    alts = alg
    dels = alts.count('*') 
    #remove insertions
    alts = _remove_indels( alts )
    #get base counts
    baseCounts = [ ( alts.upper().count(base),base ) for base in alphabet ]
    #get base frequencies
    for base_count,base in sorted( baseCounts ):
        freq = base_count*1.0/cov
        if base!=base_ref and freq >= minFreq:
            #check if alt base in both strands
            if bothStrands: 
                if not base.upper() in alts or not base.lower() in alts:
                    return
            return (base,freq) # base!=base_ref and 
  
def parse_mpileup( fnames,fastaFn,minDepth,minFreq,indels,reference,bothStrands,verbose,alphabet='ACGT' ):
    """
    """
    # open out files and write header
    header = '#coordinate\treference coverage\tref base\tref freq\talt coverage\talt base\talt freq\n'
    i = 0
    outFiles=[]
    for fn in fnames:
        i += 1
        if i==1 and reference:
            continue
        out = open( "%s.snps.cov_%s.freq_%s.bothStrands_%s.txt" % ( fn,minDepth,minFreq,bothStrands ),"w" )
        outFiles.append( out )
        out.write( header )
    
    #process mpileup
    contigs=[]
    totCov={}; totLen={}; pContig=pPos=0
    #open subprocess
    args = [ 'samtools','mpileup' ] + fnames #; print args
    proc = subprocess.Popen( args,stdout=subprocess.PIPE,bufsize=65536 )
    for line in proc.stdout:
        line      = line.strip()
        lineTuple = line.split('\t')
        #get coordinate
        refCov, refFreq = 0, 1.0
        contig,pos,baseRef = lineTuple[:3]
        samplesData = lineTuple[3:]
        #laod ref data
        if reference:
            refCov,refAlgs,refQuals = lineTuple[3:6]
            refCov = int(refCov)
            samplesData = lineTuple[6:]
            if refCov<minDepth:
                continue
            alt_allele = get_alt_allele( '',refCov,refAlgs,minFreq,alphabet,reference,bothStrands )
            if not alt_allele:
                continue
            baseRef,refFreq = alt_allele
            
        #,cov,alg,quals= #; contig,pos,base,cov,alg,quals
        for out,cov,alg,quals in zip( outFiles,samplesData[0::3],samplesData[1::3],samplesData[2::3] ):
            cov=int(cov)
            if cov<minDepth:
                continue
            # check for SNP
            alt_allele = get_alt_allele( baseRef,cov,alg,minFreq,alphabet,reference,bothStrands )
            if not alt_allele:
                continue
            # get base and freq
            base,freq = alt_allele
            lineOut='%s:%s\t%s\t%s\t%1.4f\t%s\t%s\t%1.4f\n' % (contig, pos, refCov, baseRef, refFreq, cov, base, freq)
            out.write( lineOut )
            if verbose:
                print lineOut,
  
    for out in outFiles:    #close outfile, if opened
        out.close()

def main():

    usage  = "usage: %prog [options] ref.bam *.bam" 
    parser = OptionParser( usage=usage,version="%prog 1.0" ) #allow_interspersed_args=True

    parser.add_option("-i", dest="fasta", help="fasta [required only if no reference bam]")
    parser.add_option("-d", dest="minDepth", default=5,  type=int,
                      help="minimal depth [%default]")
    parser.add_option("-f", dest="minFreq",  default=0.8, type=float,
                      help="min frequency of alternative base [%default]")
    parser.add_option("-n", dest="indels",   default=False, action="store_true", 
                      help="report indels [%default]")
    parser.add_option("-b", dest="bothStrands", default=False, action="store_true", 
                      help="only SNPs confirmed by both strand algs [%default]")
    parser.add_option("-r", dest="reference",   default=True, action="store_false", 
                      help="first bam is reference algs [%default]")
    parser.add_option("-v", dest="verbose",  default=False, action="store_true" )
  
    ( o, fnames ) = parser.parse_args()
    if o.verbose:
        print "Options: %s\nFiles: %s" % ( o,fnames )

    if not fnames:
        parser.error( "Provide at least one bam file!" )
    for fn in fnames:
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )
    
    #parse pileup
    parse_mpileup( fnames,o.fasta,o.minDepth,o.minFreq,o.indels,o.reference,o.bothStrands,o.verbose )
  
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
