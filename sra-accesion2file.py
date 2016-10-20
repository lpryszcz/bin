#!/usr/bin/env python
desc="""Retrieve SRA files for given SRA accessions.
Currently handles only experiments accessions (SRX......).

Author:
l.p.pryszcz+git@gmail.com

Barcelona, 28/09/2012
"""

import os, sys
from datetime import datetime
from optparse import OptionParser
from ftplib   import FTP

def accession2file( expacc,ftpdomain,verbose ):
    """
    """
    #connect to ftp
    ftp = FTP(ftpdomain)
    ftp.login('anonymous', '')
    #create outdir
    if not os.path.isdir( expacc ):
        os.makedirs( expacc )
    #get run accessions
    #sra/sra-instant/reads/ByExp/sra/SRX/SRX081/SRX081471/SRR305328/
    expdir = "sra/sra-instant/reads/ByExp/sra/SRX/%s/%s" % ( expacc[:6],expacc )
    dirs = ftp.nlst(expdir)
    #iterate
    if verbose:
        sys.stdout.write( "Fetching %s runs for %s...\n" % (len(dirs),expacc) )
    for dpath in dirs:
        runacc = os.path.basename( dpath )
        remotefpath = os.path.join( expdir,runacc,"%s.sra" % runacc )
        localfpath  = os.path.join( expacc,"%s.sra" % runacc )
        if not os.path.isfile( localfpath ):
            if verbose:
                sys.stdout.write( " %s > %s\n" % ( remotefpath,localfpath ) )
            #fetch file
            #continue
            lfile = open( localfpath,"w" )
            ftp.retrbinary( "RETR %s" % remotefpath, lfile.write )
            lfile.close()
        else:
            if verbose:
                sys.stdout.write( " File exists: %s\n" % ( localfpath, ) )

def main():
    usage   = "usage: %prog [options] | gzip > taxid.CS.txt.gz"
    version = "%prog 1.0"
    parser  = OptionParser( usage=usage,version=version,description=desc ) #allow_interspersed_args=True

    parser.add_option("-v", dest="verbose",  default=True, action="store_false")
    parser.add_option("-a", dest="accessions",  default="", 
                      help="comma-separated list accessions to fetch [%default]" )
    parser.add_option("-f", dest="ftp",  default="ftp-trace.ncbi.nih.gov", 
                      help="ftp server address [%default]" )
    o,args = parser.parse_args()
    if o.verbose: 
        sys.stderr.write( "Options: %s\n" % str(o) )

    #get list of accessions
    if not o.accessions:
        parser.error("Provide SRA accessions (SRX......) seperated by commas (-p)!")
    accessions = [ x for x in o.accessions.split(',') ]

    for acc in accessions:
        #get runs sra and store
        accession2file( acc,o.ftp,o.verbose )
  
if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )