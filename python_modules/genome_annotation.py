#!/usr/bin/env python
"""This file store python modules to handle 
genome annotation formats like GTF, GFF, EMBL, GenBank.
hmmer tblout etc
"""
import gzip, os, sys

def get_handle(inhandle):
    """Return opened file object
    Consider adding bzip support.
    """
    #first check if it's file or stream # hasattr(input, 'read'): #
    if type(inhandle) is file: 
        handle = inhandle
    #then check if file
    elif type(inhandle) is str:
        if not os.path.isfile(inhandle):
            sys.exit("No such file: %s" % inhandle)
        #open as gzipped if needed
        if inhandle.endswith('.gz'):
            handle = gzip.open(inhandle)
        else:
            handle =      open(inhandle)        
    else:
        sys.exit("Cannot recognise handle: %s" % inhandle)
    return handle

def _get_formatted_seq( seq,lineLen=60 ):
    """Return formatted sequence."""
    seq = str(seq)
    return "\n".join( seq[i:i+lineLen] for i in range(0,len(seq),lineLen) )

def genome2dict(inhandle, seqformat="fasta"):
    """Return genome dictionary:
    d = { contig/chromosome: seqObject, }
    Can handle fasta, genbank, etc (also gzipped).
    Requeries Bio.SeqIO (Biopython). """
    import Bio.SeqIO as SeqIO
    genome2dict = {}
    # parse genome handle
    for r in SeqIO.parse(get_handle(inhandle), seqformat):
        genome2dict[r.id] = r.seq
    return genome2dict

def get_gc( ctg2seq,contig,start=0,end=0 ):
    """Return GC of given region.
    Note: N unaware and capital letter only!"""
    if start and end:
        seq = ctg2seq[contig][start:end]
    else:
        seq = ctg2seq[contig]
    gc = 1.0 * ( seq.count("G") + seq.count("C") ) / len(seq)
    return gc

def get_contig2size( fname ):
    """Return contig2size dictionary:
    d = { contig/chromosome: size }
    Requires Biopython installed.
    """
    import Bio.SeqIO as SeqIO
    c2s = {}
    # 
    for r in SeqIO.parse( open( fname ),"fasta" ):
        contig      = r.id
        c2s[contig] = len( r.seq )
    return c2s

def get_contig2size_samtools( fname,verbose=1 ):
    """Return contig2size dictionary:
    d = { contig/chromosome: size }
    Requires samtools installed.
    """
    import os,sys
    c2s = {}
    # generate fasta index if not previously generated
    idxfn = fname+".fai"
    if not os.path.isfile( idxfn ):
        cmd = 'samtools faidx "%s" > "%s"' % ( fname,idxfn )
        if verbose:
            sys.stderr.write( cmd+"\n" )
        os.system( cmd )
    # parse fasta idex file
    for l in open( idxfn ):
        contig,size = l.split('\t')[:2]
        c2s[contig] = int( size )
    return c2s

def get_contig2readcount( fname ):
    """Return contig2size dictionary:
    d = { contig/chromosome: size }
    Requires samtools installed.
    """
    import os,sys
    c2s = {}
    # generate fasta index if not previously generated
    idxfn = fname+".idxstats"
    if not os.path.isfile( idxfn ):
        cmd = "samtools idxstats %s > %s" % ( fname,idxfn )
        sys.stderr.write( cmd+"\n" )
        os.system( cmd )
    # parse fasta idex file
    for l in open( idxfn ):
        contig,size,algs = l.split('\t')[:3]
        c2s[contig] = int( algs )
    return c2s

def get_contig2coverage( fname ):
    """Return contig2coverage dictionary:
    d = { contig/chromosome: ( readcount,size ) }
    Requires samtools installed.
    """
    import os,sys
    c2s = {}
    # generate fasta index if not previously generated
    idxfn = fname+".idxstats"
    if not os.path.isfile( idxfn ):
        cmd = "samtools idxstats %s > %s" % ( fname,idxfn )
        sys.stderr.write( cmd+"\n" )
        os.system( cmd )
    # parse fasta idex file
    for l in open( idxfn ):
        contig,size,algs = l.split('\t')[:3]
        c2s[contig] = ( int( algs ),int(size) )
    return c2s

def load_gtf(inhandle, partial=True):
    """Parse gtf and return gff info as:
    gene2position = { transcript_id: [ contig,CDSs_list,strand,function ]
    contig2gene= { contig: [ (start,end,feature,id) ]
    Note, gene2position stores CDS coordinates, while contig2gene gene boundaries!
    IMPLEMENT dealing with non-coding genes! (rRNA,tRNA,ncRNA)
    """
    import urllib
    # parse gff and take only best match for each query
    contig2gene   = {}
    gene2position = {}
    gene2boundaries = {}
    for line in get_handle(inhandle):
        line=line.strip() 
        if line.startswith('#') or not line: 
            continue
      
        line_data = line.split('\t')
        if len(line_data)<8: 
            continue # skip incorrect lines
    
        contig,source,feature,start,end,score,strand,frame,comments=line_data
    
        coordinates=[ int(start),int(end) ]
        coordinates.sort()
      
        description={}
        for atr_value in comments.split(';'):
            atr_value = atr_value.strip()
            if not atr_value:
                continue
            atr   = atr_value.split()[0]
            value = " ".join( atr_value.split()[1:] ).strip('"')
            #value = value.strip('"')
            description[atr]=value
      
        function = '' # description['Note']
        if 'note' in description:
            function=description['note']
        elif 'pfam' in description:
            function=description['pfam']
        if '%' in function:
            function = urllib.unquote(function)

        if not 'transcript_id' in description:
            #sys.stderr.write("load_gtf: Error: No transcript_id in comments: %s\n"%comments)
            continue
        id       = description['transcript_id']
        
        if feature in ('start_codon','stop_codon'):
            if id not in gene2boundaries:
                gene2boundaries[id] = []
            gene2boundaries[id] += coordinates
        elif feature == "CDS": 
            if not id in gene2position:    # add cds info to gene2position
                gene2position[id]=[ contig, [], strand, function, [] ]
            #add cds info
            gene2position[id][1].append(coordinates)
            #add frame info
            if frame.isdigit():
                gene2position[id][4].append(int(frame))
  
    # add contig2gene 
    for id in gene2position:
        contig,CDSs,strand,function,frames = gene2position[id]
        CDSs.sort()
    
        if contig not in contig2gene:
            contig2gene[contig]=[]
        
        # get gene start and stop
        ##FOR SOME REASON full transcripts are not working! likely related to
        ##start stop codons
        #print gene2boundaries[id]
        if id in gene2boundaries and len(gene2boundaries[id]) == 4: # provided correctly by gtf
            start, end = min(gene2boundaries[id]),max(gene2boundaries[id])
        elif partial:
            # take first and last CDS
            start, end = CDSs[0][0], CDSs[-1][1]
        else:
            continue
        feature = 'gene'
        gffData = ( start,end,feature,id )
        contig2gene[contig].append( gffData )
    
    # sort contig2gene
    for contig in contig2gene:
        contig2gene[contig].sort()
    return gene2position, contig2gene

def load_gff(input):
    """Parse gff and return gff info as:
    gene2position = { transcript_id: [ contig,CDSs_list,strand,function ]
    contig2gene= { contig: [ (start,end,feature,id) ]
    Note, gene2position stores CDS coordinates, while contig2gene gene boundaries!
    IMPLEMENT dealing with non-coding genes! (rRNA,tRNA,ncRNA)
    """
    import urllib
    # parse gff and take only best match for each query
    contig2gene   = {}
    gene2position = {}
    gene2boundaries = {}
    for line in get_handle(input):
        line=line.strip() 
        if line.startswith('#') or not line: 
            continue
      
        line_data = line.split('\t')
        if len(line_data)<8: 
            continue # skip incorrect lines
    
        contig,source,feature,start,end,score,strand,frame,comments=line_data
        start,end = int(start),int(end)
        coordinates=[ start,end ]
        coordinates.sort()
      
        description={}
        for atr_value in comments.split(';'):
            atr_value = atr_value.strip()
            if not atr_value:
                continue
            atr,value=atr_value.split("=") #transId=comments.split(';Name=')[-1].split(';')[0]
            value = value.strip('"')
            description[atr]=value
      
        function = '' # description['Note']
        if 'note' in description:
            function=urllib.unquote( description['note'] )
        elif 'Note' in description:
            function=urllib.unquote( description['Note'] ) #compatible with SGD
        if   "ID" in description:
            id = description['ID']
        elif "Parent" in description:
            id = description['Parent']
        else:
            sys.stderr.write( "warning: load_gtf: No ID/Parent in %s\n" % str(line_data) )

        if   feature == "gene":
            if contig not in contig2gene:
                contig2gene[contig]=[]
            gffData = ( coordinates[0],coordinates[1],feature,id )
            contig2gene[contig].append( gffData )
        elif feature == "CDS":
            if not id in gene2position:    # add cds info to gene2position
                gene2position[ id ]=[ contig,[],strand,function ]
            gene2position[ id ][1].append( coordinates ) 
  
    # sort contig2gene
    for contig in contig2gene:
        contig2gene[contig].sort()
    return gene2position, contig2gene

def load_counts_gff( fpath ):
    """Parse counts.gff and counts info for every CDS and *RNA:
    gene2position = { transcript_id: [ contig,CDSs_list,strand,function ]
    contig2gene= { contig: [ (start,end,feature,id) ]
    Note, gene2position stores CDS coordinates, while contig2gene gene boundaries!
    IMPLEMENT dealing with non-coding genes! (rRNA,tRNA,ncRNA)
    """
    import urllib
    # parse gff and take only best match for each query
    gene2count      = {}
    gene2details    = {}
    for line in open(fpath):
        line=line.strip() 
        if line.startswith('#') or not line: 
            continue
      
        line_data = line.split('\t')
        if len(line_data)<8: 
            continue # skip incorrect lines
    
        contig,source,feature,start,end,score,strand,frame,comments=line_data[:9]
        count = line_data[9]
        if not count.isdigit():
            sys.stderr.write( "load_counts_gff: Warning: Wrong read count: %s\n" % ", ".join( line_data ) )
            continue

        count       = int(count)
        start,end   = int(start),int(end)
        coordinates =[ start,end ]
        coordinates.sort()
      
        description={}
        for atr_value in comments.split(';'):
            atr_value = atr_value.strip()
            if not atr_value:
                continue
            atr,value=atr_value.split("=") #transId=comments.split(';Name=')[-1].split(';')[0]
            value = value.strip('"')
            description[atr]=value
      
        function = '' # description['Note']
        if 'note' in description:
            function=urllib.unquote( description['note'] )
        elif 'Note' in description:
            function=urllib.unquote( description['Note'] ) #compatible with SGD
        if   "ID" in description:
            id = description['ID']
        elif "Parent" in description:
            id = description['Parent']
        else:
            sys.stderr.write( "warning: load_gtf: No ID/Parent in %s\n" % str(line_data) )

        if feature == "CDS" or feature.endswith("RNA"):
            if id not in gene2details:
                gene2details[id] = [ contig,start,end,strand,feature,function ]
                gene2count[id]   = count
            else:
                gene2count[id]  += count
                if start < gene2details[id][1]:
                    gene2details[id][1] = start
                if end   < gene2details[id][2]:
                    gene2details[id][2] = end

    return gene2details,gene2count
    
def coding_snp_info( contig,geneid,CDSs,strand,ref,alt,pos ):
    """Return info about SNP:
  1. SNP type: exonic, intronic, 5'-UTR or 3'-UTR. If exonic check for:
  2. amino change type: None, synonymous, nonsynonymous, nonsense or read-through. If one of latter 3, additionally give:
  3. codon number
  4. position in codo
  5. reference codon
  6. reference amino acid
  7. SNP codon
  8. SNP amino acid
  
    """
    import Bio.Seq as Seq
    import sys
    cds = Seq.Seq('')
    cdsPos = -1
    snpType = None
    i = 0
    for s,e in CDSs:
        i += 1
        # get base position in cds
        if s<=pos<=e: 
            cdsPos = len( cds )+pos-s
            snpType = "exonic" # set exonic type
        elif pos>e:
            intronLen=0
            if len(CDSs)>i:
                intronLen = CDSs[i][0]-e
            snpType = "intronic +%s / %s" % ( pos-e,intronLen ) # set intronic - will be overwritten later if only one exon
      
        cds += contig[s-1:e]
  
    # define if UTR
    if   pos>CDSs[-1][1]:
        if strand == '-': snpType = "5'-utr"
        else:             snpType = "3'-utr"
    elif pos<CDSs[0][0]:
        if strand == '-': snpType = "3'-utr"
        else:             snpType = "5'-utr"  
  
    # reverse complement if needed
    if strand == '-':
        cds    = cds.reverse_complement()
        cdsPos = len(cds) - cdsPos - 1              # get position from other end
        ref    = str( Seq.Seq( ref ).complement() ) # get complement of ref base as well
        alt    = str( Seq.Seq( alt ).complement() ) # get complement of SNP base as well
    # check exonic or UTR SNPs
    if snpType == "exonic": # cdsPos > -1:
        # get codon & position of SNP in that codon
        iCodon    = cdsPos // 3 # get codon
        # get codon sequence
        codon     = cds[iCodon*3:iCodon*3+3] 
        codonPos  = cdsPos % 3  # position of SNP in codon
        # define alternative codon
        refCodon  = codon[:codonPos]+ref+codon[codonPos+1:] 
        altCodon  = codon[:codonPos]+alt+codon[codonPos+1:]
        # check for same codon error!
        if str(refCodon)==str(altCodon):
            sys.stderr.write( "Warning: The same codons ( %s & %s ) after SNP introduction.\n%s\n%s\n%s\n" % ( codon,altCodon,str(cds),cds[iCodon*3-6:iCodon*3+9],cds[iCodon*3-6:iCodon*3+9].translate() ) )
        # check if aa change or not
        aaType ="non-synonymous" 
        AA     = str(refCodon.translate())
        altAA  = str(altCodon.translate())
        if AA != altAA:
            if   altAA=='*':
                aaType="nonsense"
            elif AA   =='*':
                aaType="read-through"
        else:
            aaType="synonymous"
        # add info to line
        line = "%s\t%s\t%s\t%s/%s\t%s\t%s\t%s\t%s\t%s" % ( snpType,geneid,aaType,iCodon+1,len(cds)/3,codonPos,refCodon,refCodon.translate(),altCodon,altCodon.translate() )
    elif snpType:
        line = "%s\t%s\t\t\t\t\t\t\t" % (snpType,geneid)
    else:
        line = "%s\t\t\t\t\t\t\t\t" % (snpType,)
    
    # print line,# cdsPos,len(cds)
    return line 

def parse_gtf( gtf,idname='transcript_id',feature="CDS" ):
    """Parse genome annotation (GTF) and return 2 dictionaries:
    1. contig2coding={ contig: [ ( s,e,strand,transcript_id ), ] } 
    2. id2position  ={ transcript_id: [ (s,e,strand,score,frame),...  ] }
    List for each contig and transcript are sorted by start position.
    """
    import sys
    ctg2cds = {}
    id2cds  = {}
    for i,line in enumerate( open( gtf ) ):
        line = line[:-1]
        if line.startswith('#') or not line:
            continue
    
        try:  
            contig,pred,f,s,e,score,strand,frame,comment = line.split('\t')
        except:
            sys.stderr.write( "Warning: Wrong line %s: %s\n" % ( i,str(line.split('\t')) ) )
            continue
            
        f = f.upper()
        if f != feature:
            continue 
        
        #s,e,frame = int(s),int(e),float(frame)
        s,e = int(s),int(e)
        if not contig in ctg2cds:
            ctg2cds[contig] = []
        
        # get id or put fake id if no comment (ie gff2 from exonerate)
        if comment:
            k2v = {}
            for c in comment.split(';'):
                if not c or len(c.split())<2:
                    continue
                k,v = c.split()[:2]
                k,v = k.strip(), v.strip('"')
                k2v[k] = v

            # store info
            id        = k2v[ idname ]
        else:
            id = "CDS%7i" % (i+1,)
            id = id.replace(" ","0")
        #matchData = ( s,e,strand,id )
        matchData = ( s,e,id,strand,score )
        ctg2cds[contig].append( matchData )
        
        if not id in id2cds:
            id2cds[id] = []

        cdsData   = ( s,e,strand,score,frame )
        #cdsData   = ( s,e,id,strand,score )
        id2cds[id].append( cdsData )

    # sort cds for each contig/chromosome
    for ctg in ctg2cds:
        ctg2cds[ctg].sort()
    # sort cds for each transcript
    for id in id2cds:
        id2cds[id].sort()
    #print ctg2cds.keys()
    return ctg2cds,id2cds


def load_sgd_gff( fpath,test=0 ):
    """Parse sgd gff and return tuple of 3 elements:
    1. contig2coding = { chrI: [ (start,end,id,strand,score),... ] }
    2. number of genes
    3. transcript2exon_positions = { transcript_id: [ ( start,end,contig/chromosome,strand ),... ]
    NOT FINISHED!
    """
    transId2exon_pos={}
    contig2coding={};
    for line in open(fpath):
        line=line.strip() 
    
        # skip empty lines and comments
        if line.startswith('#') or not line: 
            continue
        line_data=line.split('\t')
        # skip incorrect lines
        if len(line_data)!=9: 
            continue
        contig,source,feature,start,end,score,strand,frame,comments=line_data
        # read only CDS entries
        if feature!='CDS': 
            continue
        start,end=int(start),int(end)
        # 
        transId=comments.split(';Name=')[-1].split(';')[0]
        
        # add contig to contig2position if not there 
        if contig not in contig2coding: contig2coding[contig]=[]
        
        # store coordinates on condig
        coordinates=( start,end,transId,strand,score )
        contig2coding[contig].append( coordinates ) # GeneId=comments.split('"')[1]
        
        # store gene seq
        exon_coordinates=( start,end,contig,strand )
        try:    transId2exon_pos[transId].append( exon_coordinates ) # add exon
        except: transId2exon_pos[transId]=[exon_coordinates]
        
    return contig2coding,transId2exon_pos 

def parse_decypher( infile,evalue=1e-05,qcov=0,tcov=0,verbose=True ):
    """Return hits from decypher output passing evalue and *cov
    criteria.
    matchesData = [ ( qlocus,qstart,qend,tlocus,tstart,tend,score,e ), ... ]
    """
    if verbose:
        sys.stderr.write( "Parsing decypher output...\n" )
    i = 0 
    matchData = []
    for l in open( infile ):
        l = l.strip()
        if not l or l.startswith("#"):
            continue
        lData = l.split('\t')
        if len(lData) != 15:
            sys.stderr.write( " parse_decypher: Error: Cannot parse line: %s\t" % ",".join(lData) )
            continue
        i+=1
        #get data
        qlocus,tlocus,rank,status,score,e,qstart,qend,qlength,tstart,tend,tlength,gaps,matches,sim = lData
        #get int and float
        qstart,qend,qlength,tstart,tend,tlength = int(qstart),int(qend),int(qlength),int(tstart),int(tend),int(tlength)
        score,e = float(score),float(e)
        #evalue filter
        if evalue and e>evalue:
            continue
        #tcov filter
        if tcov and tcov < (tend-tstart+1)*1.0 / tlength:
            continue
        #qcov filter
        if qcov and qcov < (qend-qstart+1)*1.0 / qlength:
            continue
        #store match info
        matchData.append( ( qlocus,qstart,qend,tlocus,tstart,tend,score,e ) )

    if verbose:
        sys.stderr.write( "  %s [%.2f%s] records passed filtering.\n" % (len(matchData),len(matchData)*100.0/i,'%') )        
    return matchData

def parse_blast( infile,q2len,t2len,evalue=1e-05,qcov=0,tcov=0,verbose=True ):
    """Return hits from decypher output passing evalue and *cov
    criteria.
    matchesData = [ ( qlocus,qstart,qend,tlocus,tstart,tend,score,e,qcoverage,tcoverage ), ... ]
    """
    import math
    if verbose:
        sys.stderr.write( "Parsing blast output...\n" )
    i = 0 
    matchData = []
    for l in open( infile ):
        l = l.strip()
        if not l or l.startswith("#"):
            continue
        lData = l.split('\t')
        if len(lData) != 12:
            sys.stderr.write( " parse_blast: Error: Cannot parse line: %s\t" % ",".join(lData) )
            continue
        i+=1
        #get data
        qlocus,tlocus,identity,algLen,mismatches,gaps,qstart,qend,tstart,tend,e,score = lData
        identity = float( identity )
        #get qlen & tlen
        qlength = tlength = 0
        if qlocus in q2len:
            qlength = q2len[qlocus]
        if tlocus in t2len:
            tlength = t2len[tlocus]
        #get int and float
        qstart,qend,qlength,tstart,tend,tlength = int(qstart),int(qend),int(qlength),int(tstart),int(tend),int(tlength)
        score,e = float(score),float(e)
        #evalue filter
        if evalue and e>evalue:
            continue
        #tcov filter
        qcoverage = tcoverage = 0
        if tlength:
            tcoverage = math.fabs( (tend-tstart+1)*1.0 / tlength )
        if tcov and tcov > tcoverage:
            continue
        #qcov filter
        if qlength:
            qcoverage = math.fabs( (qend-qstart+1)*1.0 / qlength )
        if qcov and qcov > qcoverage:
            continue
        #store match info
        matchData.append( ( qlocus,tlocus,identity,algLen,mismatches,gaps,qstart,qend,tstart,tend,e,score,qcoverage,tcoverage ) )

    if verbose:
        sys.stderr.write( "  %s [%.2f%s] records passed filtering.\n" % (len(matchData),len(matchData)*100.0/i,'%') )        
    return matchData

def parse_blat( fn,verbose=1,header=1,skipSelfMatches=False ):
    """Return list of matches from BLAT34 .psl output.
    """
    matchData=[]
    for l in open(fn):
        l=l.strip()
        if not l:
            continue
        if header:
            if l.startswith("-----"):
                header=0
            continue
        #parse line
        M,misM,repM,Ns,Qgaps,Qgapbases,Tgaps,Tgapbases,strand,Q,Qlen,Qs,Qe,T,Tlen,Ts,Te,blocks,blockSizes,Qstarts,Tstarts = l.split("\t")
        M,misM,Qgapbases,Qlen,Qs,Qe,Tlen,Ts,Te = int(M),int(misM),int(Qgapbases),int(Qlen),int(Qs),int(Qe),int(Tlen),int(Ts),int(Te)
        if skipSelfMatches and Q==T:
            continue
        #calculate missing vars
        Qalg = Qe-Qs
        Talg = Te-Ts
        identity = M*100.0/Qalg
        qcoverage = Qalg*1.0/Qlen
        tcoverage = Talg*1.0/Tlen
        #store info
        matchData.append( ( Q,T,identity,Qalg,misM,Qgapbases,Qs,Qe,Ts,Te,strand,qcoverage,tcoverage ) )
    return matchData
        
def reverse_complement( s ):
    """Return reverse complement of s.
    No error catching!
    Deprecated: Biopython is 2x faster! look under: /home/services/test/reverse_complement/"""
    d  = { "A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "-": "-",
           "a": "t", "t": "a", "c": "g", "g": "c", "n": "n",
          }
    return "".join( d[b] for b in s[::-1] )
                
def get_gap_boundaries( contig2fasta,isize=300 ):
    """Return dictionary. Chromosome/contig names are keys, and values are
    list of gap (and surrounding isize sequence) coordinates.
    Two consecutive gaps are merged if the distance between them is less than isize."""
    import re
    contig2gaps = {}
    #regex to match insert size before and after gap;
    ##matches next gap if following ACGT seq is shorter than insert size
    ###consider some stuff to recognise to big gaps maybe...
    gappat = re.compile("([ACGT]{0,%s}N+[ACGT]{0,%s})+" % (isize,isize) )
    #iterate contigs
    for c,s in contig2fasta.iteritems():
        contig2gaps[c] = []
        for m in gappat.finditer( str(s) ):
            contig2gaps[c].append( (m.start(),m.end()) )

    return contig2gaps
    
def get_edges( contig2fasta,isize=300 ):
    """Return first and last isize of each chromosome/contig."""
    contig2edge = {}
    for c,s in contig2fasta.iteritems():
        contig2edge[c] = ( isize,len(s)-isize )
    return contig2edge

def nucmer2list( fn,refSort=True ):
    """Return list of matches from nucmer.
    Sorted by ref in refSort=True."""
    matches = []
    header  = True
    for l in open( fn ):
        #skip head lines
        if l.startswith("========="):
            header = False
            continue
        if header:
            continue
        #remove white spaces from start and end
        l = l.strip()
        if not l:
            continue
        #read data
        lData = l.split()
        if   len(lData) == 13:
            rStart,rStop,sp,qStart,qStop,sp,rLen,qLen,sp,identity,sp,r,q = lData
        elif len(lData) == 11:
            rStart,rStop,sp,qStart,qStop,sp,rLen,qLen,sp,r,q = lData
            identity = 0
        else:
            sys.stderr.write( " Err: Cannot correctly parse line: %s\n" % str(lData) )
            continue
        match = ( r,int(rStart),int(rStop),q,int(qStart),int(qStop),float(identity) )
        #prepare dict
        #if r not in matches:
        #    matches[r]=[]
        matches.append( match )
        
    return matches

def nucmer2list2( fn,refSort=True ):
    """Return list of matches from nucmer.
    Sorted by ref in refSort=True."""
    matches = []
    header  = True
    for l in open( fn ):
        #skip head lines
        if l.startswith("========="):
            header = False
            continue
        if header:
            continue
        #remove white spaces from start and end
        l = l.strip()
        if not l:
            continue
        #read data
        rStart,rStop,sp,qStart,qStop,sp,rLen,qLen,sp,identity,sp,r,q = l.split()
        match = ( q,int(qStart),int(qStop),r,int(rStart),int(rStop),float(identity) )
        #prepare dict
        #if r not in matches:
        #    matches[r]=[]
        matches.append( match )
        
    return matches

def lastal2list( fn,refSort=True ):
    """Return list of matches from lastal tab formatted file.
    """
    matches = []
    for l in open( fn ):
        l = l.strip()
        if not l or l.startswith("#"):
            continue
        score,r,rStart,rAlnSize,strand1,seqSize1,q,qStart,qAlnSize,strand2,seqSize2,blocks = l.split()
        rStart = int(rStart)
        rStop  = rStart + int(rAlnSize)
        qStart = int(qStart)
        qStop  = qStart + int(qAlnSize)
        match = ( r,rStart,rStop,q,qStart,qStop,float(score) )

        matches.append( match )
    #sort matches by reference start
    matches.sort()
        
    return matches 

def load_transcripts_bed( bed,oneoff=False ):
    """Return dictionary of transcripts from BED.
    transcripts = { transcript: { strand: +/-, intervals: [] } }
    If oneoff, add +1 to start and end coordinate."""
    transcripts = {}
    for l in open(bed):
        l = l.strip()
        if not l or l.startswith("#"):
            continue
        ref,start,end,name,score,strand = l.split("\t")[:6]
        start,end = int(start),int(end)
        if oneoff:
            start+=1
            end  +=1
        if not name in transcripts:
            transcripts[name]={ "chromosome": ref, "strand": strand, "intervals": [] }
        transcripts[name]["intervals"].append( (start,end,score) )
    return transcripts

    