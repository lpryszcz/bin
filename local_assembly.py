#!/usr/bin/env python
"""Run local assembly pipeline.

Prerequisites:
Ray
samtools

Author:
l.p.pryszcz@gmail.com

Barcelona, 16/05/2012
"""

import os, sys
import gzip
import subprocess
from optparse import OptionParser,OptionGroup
from datetime import datetime
from genome_annotation import get_contig2readcount

def get_read_ids( bam,contigs,outdir,split,positions,gzipped,force,verbose ):
    """Get read ids of aligned reads."""
    idFnames = []
    if verbose:
        sys.stderr.write( "Getting mapped read ids...\n" )    
    ids = set()
    if not split:
        outfn = os.path.join(outdir,"selected.ids")
        if gzipped:
            outfn += ".gz"
        idFnames.append( outfn )
        #skip if not force and file exists    
        if not force and os.path.isfile( outfn ):
            return idFnames
        #open outfile
        if gzipped:
            out = gzip.open( outfn,"w" )
        else:
            out =      open( outfn,"w")
    #process chromosomes
    for c in contigs:
        if verbose:
            sys.stderr.write( " %s           \r" % c )
        #skip *
        if c=="*":
            continue
        #if split by chr
        if split:
            outfn = os.path.join(outdir,c+".ids")
            if gzipped:
                outfn += ".gz"
            idFnames.append( outfn )
            #skip if not force and file exists    
            if not force and os.path.isfile( outfn ):
                continue
            #empty set and open out file
            ids   = set()
            #open outfile
            if gzipped:
                out = gzip.open( outfn,"w" )
            else:
                out =      open( outfn,"w")
        #launch samtools subprocess
        args = [ "samtools","view","-f1","-F4",bam,c] #,"|","cut","-f1" ]
        proc = subprocess.Popen( args,stdout=subprocess.PIPE )
        #process s and store to file
        for l in proc.stdout:
            l=l.strip()
            id = l.split('\t')[0]
            #if positions instead of read names
            if positions:
                #store int (use less space than string in hash table)
                id = int(id)
            if id and id not in ids:
                ids.add( id )
                out.write( "%s\n" % id )
        #close out for ids splitted by chr
        if split:
            out.close()
    #close out & return
    if not split:
        out.close()
        return idFnames
        
    #add unaligned pairs
    if verbose:
        sys.stderr.write( "Getting unmapped pairs ids...\n" )
    outfn = os.path.join(outdir,"unmapped.ids")
    if gzipped:
        outfn += ".gz"
    idFnames.append( outfn )
    #skip if not force and file exists    
    if not force and os.path.isfile( outfn ):
        return idFnames
    ids  = set()
    #open outfile
    if gzipped:
        out = gzip.open( outfn,"w" )
    else:
        out =      open( outfn,"w")
    args = [ "samtools","view","-f12",bam ]
    proc = subprocess.Popen( args,stdout=subprocess.PIPE )
    #process s and store to file
    for l in proc.stdout:
        l=l.strip()
        id = l.split('\t')[0]
        #if positions instead of read names
        if positions:
            #store int (use less space than string in hash table)
            id = int(id)
        if id and id not in ids:
            ids.add( id )
            out.write( "%s\n" % id )
    #close unmapped file
    out.close()
    return idFnames

def fqfetch_by_id( idsFn,fqFn,pairI,addPairInfo,gzipped,force=False ):
    """Store reads for given ids.
    Add nStart and the beginning of each read name, and
    nEnd at the end.
    Note, it's not tested yet!"""
    outFn = idsFn+"_%s.fastq" % pairI
    if gzipped:
        outFn += ".gz"
    if not force and os.path.isfile( outFn ):
        return
    #get formatting
    nStart,nEnd = "@",""
    if addPairInfo:
        nEnd = "/%s" % pairI
    #open fqFn & idsFile
    if fqFn.endswith( ".gz" ):
        handle = gzip.open( fqFn  )
        idsFile= gzip.open( idsFn )
    else:
        handle =      open( fqFn  )
        idsFile=      open( idsFn )
    #get positions of interest
    ids    = set( nStart+x.split()[0]+nEnd for x in idsFile )
    #open outfile
    if gzipped:
        out = gzip.open( outFn,"w" )
    else:
        out =      open( outFn,"w" )    
    #get reads of interest
    lineno = 0
    flag   = False
    for l in handle:
        lineno += 1
        if lineno%4 == 1: 
            flag = (l.split()[0] in ids)
        if flag:
            out.write( l )
    #close output
    out.close()

def fqfetch_by_position( idsFn,fqFn,pairI,addPairInfo,gzipped,force=False ):
    """Store reads for given positions."""
    outFn = idsFn+"_%s.fastq" % pairI
    if gzipped:
        outFn += ".gz"
    if not force and os.path.isfile( outFn ):
        return
    #open fqFn & idsFile
    if fqFn.endswith( ".gz" ):
        handle = gzip.open( fqFn  )
        idsFile= gzip.open( idsFn )
    else:
        handle =      open( fqFn  )
        idsFile=      open( idsFn )
    #get positions of interest
    #note in BAM positions are 1-based, here we prefer 0-based
    pos = set( int(x)-1 for x in idsFile )
    #open outfile
    if gzipped:
        out = gzip.open( outFn,"w" )
    else:
        out =      open( outFn,"w" ) 
    #get reads of interest
    i   = 0
    for l in handle:
        if i/4 in pos:
            out.write( l )
        i += 1
    #close output        
    out.close()
    
def main():
    usage  = "usage: %prog [options] fq1 fq2"
    usage += "\nzcat fastq | bowtie -S --all -n0 -l28 bowtie/markers.fna -q - | samtools view -SbuF 4 - | %prog [options] fq1 fq"
    desc   = """BAM file has to be sorted and indexed.
If you are interested in assembly of given chromosomes/contigs,
it's recommended to split reads by chromosome/contig (-s).
If you want to assemble some marker surrounding, leave default settings.
    """
    epilog = "Note, it may require substantial amoung of RAM memory for read ids storing."
    parser = OptionParser( usage=usage,version="%prog 1.0",description=desc,epilog=epilog ) 

    parser.add_option("-i", dest="infile",  default="",
                      help="bam file [sam from stdin]")
    parser.add_option("-o", dest="outdir",  default="local_assembly",
                      help="output directory [%default]")
    parser.add_option("-s", dest="split",   default=False, action="store_true",
                      help="split by chromosome/contig [%default]")
    parser.add_option("-p", dest="pos",     default=False, action="store_true",
                      help="position rather than read name are given [%default]")
    parser.add_option("-t", dest="addPairInfo", default=False, action="store_true",
                      help="add /1 and /2 to read ids [%default]")
    parser.add_option("-f", dest="force",   default=False, action="store_true",
                      help="overwrite        [%default]")    
    parser.add_option("-z", dest="gzipped", default=False, action="store_true",
                      help="gzip all files   [%default]")
    parser.add_option("-v", dest="verbose", default=False, action="store_true" )

    rayo = OptionGroup(parser, "Assembler options")
    rayo.add_option("-k", dest="kmer",      default="31", type=int,
                      help="K-mer size       [%default]")
    parser.add_option_group(rayo)

  
    ( o, fnames ) = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\nArgs: %s\n" % ( o,fnames ) )

    #check files
    if len(fnames)!=2:
        parser.error( "Provide 2 fastq files as arguments!" )
    for fn in [ o.infile, ] + fnames:
        if not fn:
            parser.error( "Provide input file!" )
        if not os.path.isfile( fn ):
            parser.error( "No such file: %s" % fn )

    #create output directory
    if not os.path.isdir( o.outdir ):
        os.makedirs( o.outdir )
    else:
        sys.stderr.write( "WARNING: Output directory: %s exists!\n" % o.outdir )

    #get all ref chromosomes/contigs
    if o.verbose:
        sys.stderr.write( "Getting chromosomes...\n" )        
    contig2readcount = get_contig2readcount( o.infile )
    contigs = sorted( contig2readcount.keys() )

    #select ids
    idFnames = get_read_ids( o.infile,contigs,o.outdir,o.split,o.pos,o.gzipped,o.force,o.verbose )

    #fetch reads from fq files
    if o.verbose:
        sys.stderr.write( "Fetching reads...\n" )
    for idsFn in idFnames:
        i=0
        for fqFn in fnames:
            if o.verbose:
                sys.stderr.write( " for %s from %s        \r" % (idsFn,fqFn) )        
            #to deal with various fastq formats
            i+=1
            #fetch fastq by position or id
            if o.pos:
                fqfetch_by_position( idsFn,fqFn,i,o.addPairInfo,o.gzipped,o.force )
            else:
                fqfetch_by_id( idsFn,fqFn,i,o.addPairInfo,o.gzipped,o.force )
    
if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )


