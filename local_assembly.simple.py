#!/usr/bin/env python2
"""Run local assembly using subset of reads aligning to particular sequence(s).

Prerequisites:
 Ray assembler  (http://sourceforge.net/projects/denovoassembler/)
 bowtie aligner (http://bowtie-bio.sourceforge.net/ ; version 1)
 zlib [to enable -z]

Author:
l.p.pryszcz+git@gmail.com

Barcelona, 16/05/2012
"""

import os, sys
import commands
import gzip
import subprocess
from optparse import OptionParser,OptionGroup
from datetime import datetime
from Bio      import Seq

def _get_bowtie_proc( fn,ref,mismatches,seedlen,gzipped,verbose,bufsize=65536 ):
    """Return bowtie subprocess."""
    catArgs = ['cat',fn]
    bwtArgs = ['bowtie','-Sq','--quiet','-k 1','-n',str(mismatches),'-l',str(seedlen),ref,'-']
    if gzipped:
        catArgs[0] = 'z'+catArgs[0]
    if verbose:
        sys.stderr.write( "  %s\n" % " | ".join( " ".join(args) for args in (catArgs,bwtArgs) ) )
    #select ids
    catProc = subprocess.Popen( catArgs,bufsize=bufsize,stdout=subprocess.PIPE )
    bwtProc = subprocess.Popen( bwtArgs,bufsize=bufsize,stdin=catProc.stdout,stdout=subprocess.PIPE )
    return bwtProc
            
def sam2fastq( sam ):
    """Return tuple of fastq and bool of sam entry.
    Bool is true if sam entry is correctly aligned."""
    fq,sig = "",False
    #for reference look at http://samtools.sourceforge.net/samtools.shtml
    qname,flag,tname,pos,mapq,cigar,mrnm,mpos,tlen,seq,qual = sam.split("\t")[:11]
    #get binary flag
    bflag = bin(int(flag))
    #if mapped to - strand (check if bflag has adequate length first)
    if len(bflag)>6 and bflag[-5]=='1':
        #get reverse complement of seq
        seq   = str( Seq.Seq( seq ).reverse_complement() )
        qual  = qual[::-1]
    #get fastq
    fq = "@%s\n%s\n+\n%s\n" % ( qname,seq,qual )
    #get alg signal
    if tname!="*":
        sig=True
    return ( fq,sig )
    
def get_aligning_pairs( fnames,ref,outdir,mismatches,seedlen,gzipped,force,verbose ):
    """Save pairs for which at least one read has valid alg.
    """
    outfn1 = os.path.join(outdir,"reads_1.fastq")
    outfn2 = os.path.join(outdir,"reads_2.fastq")
    if gzipped:
        outfn1 += ".gz"
        outfn2 += ".gz"
    #if not force and file exists    
    if not force and os.path.isfile( outfn1 ) and os.path.isfile( outfn2 ):
        if verbose:
            sys.stderr.write( " files exists: %s & %s\n" % (outfn1,outfn2) )
        return 
    #open outfile
    if gzipped:
        out1 = gzip.open( outfn1,"w" )
        out2 = gzip.open( outfn2,"w" )        
    else:
        out1 =      open( outfn1,"w" )
        out2 =      open( outfn2,"w" )        

    #create bowtie index
    idxfn = ref + ".1.ebwt"
    if not os.path.isfile( idxfn ):
        cmd = "bowtie-build %s %s" % (ref,ref)
        if verbose:
            sys.stderr.write( " Creating index: %s...\n" % cmd )
        bwtmessage = commands.getoutput( cmd )

    if verbose:
        sys.stderr.write( " Aligning...\n" )
    #process both pe files
    ids = set()
    while fnames:
        #get pairs fnames
        fn2 = fnames.pop()
        fn1 = fnames.pop()
        #run bowtie as subprocess
        proc1 = _get_bowtie_proc( fn1,ref,mismatches,seedlen,gzipped,verbose )
        proc2 = _get_bowtie_proc( fn2,ref,mismatches,seedlen,gzipped,verbose )        
        #process both bowtie output one by one
        k=pk=l=l1=l2=0
        while True:
            sam1 = proc1.stdout.readline()
            sam2 = proc2.stdout.readline()
            #break if finished
            if not sam1 or not sam2:
                break
            #skip sam header info
            if sam1.startswith("@"):
                continue
            fq1,sig1 = sam2fastq(sam1)
            fq2,sig2 = sam2fastq(sam2)            
            #store fq
            if sig1 or sig2:
                out1.write( fq1 )
                out2.write( fq2 )
                l  += 1
            if sig1:
                l1 += 1
            if sig2:
                l2 += 1
            k += 1
            #print progress info
            if k>pk:
                pk+=10**5
                if verbose:
                    sys.stderr.write( "   %9i [matched: %6.2f%s]    \r" % (k,l*100.0/k,'%') )
        #print summary for given reads
        if verbose:
            sys.stderr.write( "#Processed %s pairs.\n %s aligned [%.2f%s]\n %s first from pair aligned [%.2f%s]\n %s second from pair aligned [%.2f%s]\n" % ( k,l,l*100.0/k,'%',l1,l1*100.0/k,'%',l2,l2*100.0/k,'%' ) )
           
def main():
    usage  = "usage: %prog [options] fqA1 fqA2 [fqB1 fqB2]"
    desc   = """Retrieve reads (and their pairs) reads that align on given reference sequence.
Subsequently, assembly is created with given subset of reads. This (often) nicely recover
regions that are difficult to assemble.
Reference sequence can be: some marker(s), mitochondrion genome or even chromosome(s).
It takes ~20 min for 11 mln reads (2x46bp) and mitochondrion reference of 30kb.
    """
    epilog = ""
    parser = OptionParser( usage=usage,version="%prog 1.0",description=desc,epilog=epilog ) 

    parser.add_option("-o", dest="outdir",  default="simple_assembly",
                      help="output directory       [%default]")
    parser.add_option("-f", dest="force",   default=False, action="store_true",
                      help="overwrite              [%default]")    
    parser.add_option("-z", dest="gzipped", default=False, action="store_true",
                      help="gzip all files         [%default]")
    parser.add_option("-v", dest="verbose", default=False, action="store_true" )

    bwto = OptionGroup(parser, "Aligner options")
    bwto.add_option("-l", dest="seedlen",   default="28", type=int,
                      help="length of seed         [%default]")
    bwto.add_option("-n", dest="mismatches",default="2", type=int,
                      help="max mismatches in seed [%default]")
    bwto.add_option("-r", dest="ref",      default="",
                      help="reference fasta        [mandatory]")
    parser.add_option_group(bwto)

    rayo = OptionGroup(parser, "Assembler options")
    rayo.add_option("-a", dest="runAssembly", default=True, action="store_false",
                    help="run assembly             [%default]")
    rayo.add_option("-k", dest="kmer",      default="31", type=int,
                      help="K-mer size             [%default]")
    rayo.add_option("-p", dest="peakCov",   default="4000", type=int,
                      help="peak assembly coverage [%default]")
    rayo.add_option("-m", dest="minCov",    default="10", type=int,
                      help="min. assembly coverage [%default]")
    parser.add_option_group(rayo)
  
    ( o, fnames ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nFastQ files: %s\n" % ( o,fnames ) )

    #check files
    if len(fnames) % 2 != 0:
        parser.error( "Provide even number of fastq files as arguments!" )
    for fn in fnames + [ o.ref, ]:
        if not fn:
            parser.error( "Provide input file!" )
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )

    #create output directory
    if not os.path.isdir( o.outdir ):
        os.makedirs( o.outdir )
    else:
        sys.stderr.write( "WARNING: Output directory: %s exists!\n" % o.outdir )

    #get read ids
    if o.verbose:
        sys.stderr.write( "Fetching reads...\n" )        
    get_aligning_pairs( fnames,o.ref,o.outdir,o.mismatches,o.seedlen,o.gzipped,o.force,o.verbose )

    if not o.runAssembly:
        return
    #run Ray
    if o.verbose:
        sys.stderr.write( "Running assembly...\n" )
    rayDir = os.path.join( o.outdir,"ray.k%s" % o.kmer )        
    if os.path.isdir( rayDir ):
        if o.force:
            os.system( "rm -r %s" % rayDir )
        else:
            return
            
    #define fq file name base
    fqfname = os.path.join( o.outdir,"reads_?.fastq" )
    if o.gzipped:
        fqfname += ".gz"
    cmd = "Ray -k %s -p %s -o %s -peakCoverage %s -minimumCoverage %s > %s.log" % ( o.kmer,fqfname,rayDir,o.peakCov,o.minCov,rayDir )
    if o.verbose:
        sys.stderr.write( " %s\n" % cmd )
    os.system( cmd )

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )


