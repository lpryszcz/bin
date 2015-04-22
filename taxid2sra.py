#!/usr/bin/env python
desc="""Fetch all entries from SRA for given taxid.
Save the biggest run per each SAMPLE (SRS) from given date. Paired first, if any. 

Note, it run fastq-dump in background. Make sure you have enough free cores;)

DEPENDENCIES:
Biopython
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 2/10/2012
"""
changelog="""
1.1:
- argparse
- 
"""

import argparse, os, sys
from datetime import datetime
#from optparse import OptionParser
from ftplib   import FTP
from Bio      import Entrez
import xml.etree.ElementTree as ET

def srr2info( srr ):
    """Return info for SRR entry
     - experiment id
     - submission id
     - project id
     - biosample id
     - run date
     - bases
     - insert size
     - insert std
     - reads orientation
    """
    
    '''
for child in root[0]: print child.tag, child.attrib
EXPERIMENT {'center_name': 'BI', 'alias': '74116.WR23613.Solexa-42619.62C7UAAXX100916.P', 'accession': 'SRX026545'}
SUBMISSION {'submission_date': '2009-06-01T02:01:25Z', 'lab_name': 'Genome Sequencing', 'submission_comment': 'Produced by user cristyn on Sun May 31 22:01:25 EDT 2009', 'alias': 'BI.Streptococcus_pyogenes_Pathogenomics', 'center_name': 'BI', 'accession': 'SRA008647'}
STUDY {'center_name': 'BI', 'alias': 'Fusarium_oxysporum_Diversity_RNA_Sequencing_multi_isolate', 'accession': 'SRP002351'}
SAMPLE {'center_name': 'BI', 'alias': '74336.0', 'accession': 'SRS190364'}
RUN_SET {}

    root[0][0].keys()
    ['center_name', 'alias', 'accession']
    '''
    #search NCBI
    result = Entrez.read( Entrez.esearch(db="sra",term=srr ) )
    if   not result['IdList']:
        sys.stderr.write( "  Entrez Error: No results for %s\n" % srr )
        return
    elif len(result['IdList'])>1:
        sys.stderr.write( "  Entrez Warning: Multiple hits for %s: %s\n" % (srr,",".join(result['IdList'])) )
    #fetch info from NCBI
    xml = Entrez.efetch( db="sra",id=result['IdList'][0] ).read()
    root   = ET.fromstring(xml)#; print xml
    #get experiment
    EXPERIMENT = root[0].find("EXPERIMENT")
    srx = EXPERIMENT.attrib['accession']    
    #get submission
    s   = root[0].find("SUBMISSION")
    sra = s.attrib['accession']
    #get accession
    s   = root[0].find("STUDY")
    srp = s.attrib['accession']
    #get accession
    s   = root[0].find("SAMPLE")
    srs = s.attrib['accession']
    
    ''' <Element 'RUN_SET' at 0x31a1590>
RUN{'run_date': '2010-09-16T04:00:00Z',
    'center_name': 'BI',
    'total_bases':  '3035446384',
    'run_center': 'BI',
    'accession': 'SRR190806',
    'total_spots': '19970042',
    'cluster_name': 'public',
    'alias': 'BI.PE.100916_SL-XDX_00019_FC62C7UAAXX.7.srf',
    'instrument_name': 'SL-XDX',
    'published': '2011-04-28 14:48:11',
    'static_data_available': '1',
    'is_public': 'true',
    'load_done': 'true',
    'size': '1839421668'}
    '''
    #run data
    s     = root[0].find('RUN_SET') #it's within RUN_SET
    date  = s[0].attrib['run_date'] 
    bases = s[0].attrib['total_bases']

    #LIBRARY_LAYOUT - maybe try to simplify it
    isize=istdv=orient = 0
    DESIGN = EXPERIMENT.find("DESIGN") # [2][2][4][0].attrib#; print layout
    LIBRARY_DESCRIPTOR = DESIGN.find("LIBRARY_DESCRIPTOR")
    LIBRARY_LAYOUT = LIBRARY_DESCRIPTOR.find("LIBRARY_LAYOUT")
    PAIRED = LIBRARY_LAYOUT.find("PAIRED")
    if PAIRED is not None:
        layout = PAIRED.attrib
        isize  = layout['NOMINAL_LENGTH'] # NOMINAL_LENGTH="476"
        orient = layout['ORIENTATION'] #ORIENTATION="5\'3\'-3\'5\'
        istdv  = layout['NOMINAL_SDEV'] ##PAIRED NOMINAL_SDEV="149.286"
    return ( srx,sra,srp,srs,date,bases,isize,istdv,orient )

def xml2data( child,taxid2srs,verbose ):
    """ """
    #get experiment
    EXPERIMENT = child.find("EXPERIMENT")
    srx = EXPERIMENT.attrib['accession']
    #get submission
    s   = child.find("SUBMISSION")
    sra = s.attrib['accession']
    #get accession
    s   = child.find("STUDY")
    srp = s.attrib['accession']
    #get accession
    for SAMPLE in child.findall("SAMPLE"):
        #if SAMPLE.attrib['accession']!=
        srs = SAMPLE.attrib['accession']
        #get taxid
        SAMPLE_NAME = SAMPLE.find("SAMPLE_NAME")
        TAXON_ID    = SAMPLE_NAME.find("TAXON_ID")
        taxid       = int(TAXON_ID.text)
        SCIENTIFIC_NAME = SAMPLE_NAME.find("SCIENTIFIC_NAME")
        #malformed xml?
        if SCIENTIFIC_NAME is None:
            return taxid2srs
        strain          = SCIENTIFIC_NAME.text
        #get strain tag - this may cause problems with non-ENA accessions!
        SAMPLE_ATTRIBUTES = SAMPLE.find("SAMPLE_ATTRIBUTES")
        if SAMPLE_ATTRIBUTES is None:
            continue
        for SAMPLE_ATTRIBUTE in SAMPLE_ATTRIBUTES.findall("SAMPLE_ATTRIBUTE"):
            #print SAMPLE_ATTRIBUTE.find("TAG").text
            if SAMPLE_ATTRIBUTE.find("TAG").text == "strain":
                #print SAMPLE_ATTRIBUTE.find("VALUE")
                strain += " %s" % SAMPLE_ATTRIBUTE.find("VALUE").text
                break
        if strain!="unidentified organism":
            break
        
    #LIBRARY_LAYOUT - maybe try to simplify it
    isize=istdv=orient = 0
    DESIGN = EXPERIMENT.find("DESIGN") # [2][2][4][0].attrib#; print layout
    LIBRARY_DESCRIPTOR = DESIGN.find("LIBRARY_DESCRIPTOR")
    LIBRARY_LAYOUT = LIBRARY_DESCRIPTOR.find("LIBRARY_LAYOUT")
    PAIRED = LIBRARY_LAYOUT.find("PAIRED")
    if PAIRED is not None:
        layout = PAIRED.attrib
        isize  = 0
        if 'NOMINAL_LENGTH' in layout:
            isize = float( layout['NOMINAL_LENGTH'] ) # NOMINAL_LENGTH="476"
        elif not verbose:
            sys.stderr.write( "  %s: Paired run and no or zero insert. Skipped!\n" % srx ) 
        istdv  = 0
        if 'NOMINAL_SDEV' in layout:
            istdv  = float( layout['NOMINAL_SDEV'] ) ##PAIRED NOMINAL_SDEV="149.286"
        elif verbose:
            sys.stderr.write( "  %s: Paired run and no NOMINAL_SDEV\n" % srx )
        orient = ""
        if 'ORIENTATION' in layout:
            orient = layout['ORIENTATION'] #ORIENTATION="5\'3\'-3\'5\'
        elif verbose:
            sys.stderr.write( "  %s: Paired run and no orientation\n" % srx )
        
    #run data
    runs = []
    RUN_SET = child.find('RUN_SET') #it's within RUN_SET
    for RUN in RUN_SET.findall("RUN"):
        srr   = RUN.attrib['accession']
        date  = ""
        if 'run_date' in RUN.attrib:
            date  = RUN.attrib['run_date']
        elif verbose:
            sys.stderr.write( "  %s (%s): No run date!\n" % (srx,srr) )
        bases = 0
        if 'total_bases' in RUN.attrib:
            bases = int( RUN.attrib['total_bases'] )
        elif verbose:
            sys.stderr.write( "  %s (%s): No total bases!\n" % (srx,srr) )            
        runs.append( (srr,bases,date) )
        
    #store data
    childdata = ( strain,taxid,srx,srp,isize,istdv,orient,runs )
    if verbose:
        sys.stderr.write( "  %s: %s: %s\n" % (taxid,srs,str(childdata)) )
    if not taxid in taxid2srs:
        taxid2srs[taxid] = {}
    if not srs   in taxid2srs[taxid]:
        taxid2srs[taxid][srs] = []
    taxid2srs[taxid][srs].append( childdata )
    return taxid2srs
    
def taxid2runs( outfn,taxid,verbose,db="sra",retmode="xml",retmax=10**6 ):
    """Return info from SRA for given taxid. """
    taxid2srs = {}
    #search NCBI
    term   = "txid%s[organism] AND dna_data[filter] AND sra_public[filter] AND type_genome[filter]" % taxid
    if verbose:
        sys.stderr.write( "Query: %s\n" % term )
    result = Entrez.read( Entrez.esearch( db=db,term=term,retmax=retmax ) )
    ids    = result['IdList']
    if   not ids:
        sys.stderr.write( "  Entrez Error: No results for %s\n" % taxid )
        return
    if verbose:
        sys.stderr.write( "Downloading %s entries from NCBI %s database...\n" % ( len(ids),db ) )
    #post NCBI query
    '''search_handle     = Entrez.epost( db,id=",".join( ids ) )
    search_results    = Entrez.read( search_handle )
    webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"] 
    #fetch info from NCBI
    xml    = Entrez.efetch( db=db,retmode=retmode,retmax=retmax,webenv=webenv,query_key=query_key ).read()    
    root   = ET.fromstring( xml )
    for child in root:
        #update dict
        srs2exp = xml2data( child,srs2exp,verbose )'''
    for id in ids:
        xml     = Entrez.efetch( db=db,retmode=retmode,id=id ).read()#; print xml
        root    = ET.fromstring( xml )
        child   = root[0]
        taxid2srs = xml2data( child,taxid2srs,verbose )

    #print output
    out = open( outfn,"w" )
    header = "#Strain\tTaxid\tSample\tExperiment\tProject\tInsert size\tOrientation\tRun\tBases\tDate\n"
    out.write( header )
    sys.stderr.write( "Saving SRA info to: %s\n" % outfn )
    for taxid in taxid2srs:
        for srs in taxid2srs[taxid]:
            for strain,taxid,srx,srp,isize,istdv,orient,runs in taxid2srs[taxid][srs]:
                for srr,bases,date in runs:
                    line = u"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (strain,taxid,srs,srx,srp,isize,orient,srr,bases,date)
                    out.write( line.encode('ascii', 'xmlcharrefreplace') )
    out.close()

    return taxid2srs
    
def srr2fastq( odir,srr,srs,srx,ncbitimestamp,isize,ftpdomain,verbose ):
    """Fetch file."""
    #connect to ftp
    ftp = FTP(ftpdomain)
    ftp.login('anonymous', '')
    #create outdir
    outdir = os.path.join( odir,srs,srx )
    if not os.path.isdir( outdir ):
        os.makedirs( outdir )
    #get nice date timestamp
    date   = datetime.strptime( ncbitimestamp,"%Y-%m-%dT%H:%M:%SZ" )
    tstamp = date.strftime("%Y%m%d")
    #ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR540/SRR540367
    remotefpath = "sra/sra-instant/reads/ByRun/sra/SRR/%s/%s/%s.sra" % ( srr[:6],srr,srr )
    localfpath  = os.path.join( outdir,"%s_%s.sra" % ( tstamp,srr) )
    if not os.path.isfile( localfpath ):
        if verbose:
            sys.stdout.write( " %s > %s\n" % ( remotefpath,localfpath ) )
        #fetch file
        lfile = open( localfpath,"w" )
        ftp.retrbinary( "RETR %s" % remotefpath, lfile.write )
        lfile.close()
    else:
        if verbose:
            sys.stdout.write( " File exists: %s\n" % ( localfpath, ) )
    #get paired if insert size fastq
    if isize:
        cmd = "fastq-dump --gzip --split-3 -O %s %s 2>&1 > %s.log &" % ( outdir,localfpath,localfpath )
        fqfn1 = "%s_1.fastq.gz" % localfpath[:-4]
        fqfn2 = "%s_2.fastq.gz" % localfpath[:-4]
        if not os.path.isfile(fqfn1) or not os.path.isfile(fqfn2):
            if verbose:
                sys.stderr.write( "  %s\n" % cmd )
                os.system( cmd )
    else:
        cmd = "fastq-dump --gzip -O %s %s 2>&1 > %s.log &" % ( outdir,localfpath,localfpath )
        fqfn1 = "%s.fastq.gz" % localfpath[:-4]
        if not os.path.isfile(fqfn1): # or not os.path.isfile(fqfn2):
            if verbose:
                sys.stderr.write( "  %s\n" % cmd )
                os.system( cmd )        
    
def get_runs( taxid2srs,ftpdomain,orientth,maxisize,paired,minbases,verbose ):
    """Select the best run for each uniq taxid-srs-date combination
    """
    if verbose:
        sys.stderr.write( "Fetching best run for each uniq taxid-srs-date combination...\n" )
    #select the best run for each uniq taxid-srs-date combination
    for taxid in taxid2srs:
        for srs in taxid2srs[taxid]:
            date2runs={}
            for strain,taxid,srx,srp,isize,istdv,orient,runs in taxid2srs[taxid][srs]:
                #check if paired
                if paired:
                    if not isize:
                        continue
                    #skip if wrong orientation
                    if orientth and orientth!=orient:
                        continue   
                #skip big insert size or not paired
                if maxisize:
                    if isize>maxisize:
                        continue
                #add runs passed filtering   
                for srr,bases,date in runs:
                    #skip if too small yield
                    if bases < minbases*10**6:
                        continue
                    if date not in date2runs:
                        date2runs[date]=[]
                    date2runs[date].append( (srr,srx,srp,isize,bases) )
            #process best run for each uniq taxid-srs-date combination
            for date in date2runs:
                #
                fltruns = filter( lambda x: x[3]!=0, date2runs[date] )
                if not fltruns:
                    fltruns = date2runs[date]
                #sort by size
                bestrun = sorted( fltruns,key=lambda x: x[-1],reverse=True )[0]
                #print bestrun,date2runs[date]
                srr,srx,srp,isize,bases = bestrun
                #fetch
                odir = "taxid%s" % taxid
                srr2fastq( odir,srr,srs,srx,date,isize,ftpdomain,verbose )
                
def main():
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.1') 
    parser.add_argument("-d", "--download",    default=False, action="store_true", 
                        help="download SRA files")        
    parser.add_argument("-t", "--taxid",       type=int, required=True,
                        help="taxid of interest " )
    parser.add_argument("-f", dest="ftp",      default="ftp-trace.ncbi.nih.gov", 
                        help="ftp server address [%(default)s]" )
    parser.add_argument("-e", "--email",       default="lpryszcz@crg.es", type=str,
                        help="email address      [%(default)s]" )
    parser.add_argument("-o", dest="orient",   default="5'3'-3'5'",
                        help="orientation        [%(default)s]" )
    parser.add_argument("-m", dest="maxisize", default=1000, type=int,
                        help="max allowed insert [%(default)s]" )
    parser.add_argument("-b", dest="minbases", default=600, type=int,
                        help="min Mbases in run  [%(default)s Mbases -> 10x for 60Mb genome]" )
    parser.add_argument("-p", "--paired",      default=False, action="store_true",
                        help="fetch only paired runs" )
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write( "Options: %s\n" % str(o) )

    Entrez.email = o.email

    #get all runs for taxid
    outfn = "sra.tsv"
    taxid2srs = taxid2runs( outfn,o.taxid,o.verbose )

    if o.download:
        #fetch best srr
        get_runs( taxid2srs,o.ftp,o.orient,o.maxisize,o.paired,o.minbases,o.verbose )
  
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )