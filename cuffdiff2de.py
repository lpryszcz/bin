#!/usr/bin/env python
desc="""Report differentially expressed genes and their annotation.

v1.1:
- added expression for all controls, as this cause misunderstandings sometimes
"""
epilog="""AUTHOR:
l.p.pryszcz+git@gmail.com

Barcelona, 17/07/2012
"""

import argparse, os, sys, urllib
from datetime import datetime
from pfam     import load_pfam_tblout

def load_gff_annotation(handle):
    """Parse sgd gff and return tuple of 3 elements:
    """
    id2ann = {}
    for line in handle:
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
        if feature!='gene': 
            continue
        id = comments.split(';Name=')[-1].split(';')[0]
        note = comments.split(';Note=')[-1].split(';')[0]
        id2ann[id] = (urllib.unquote_plus(note),)
        
    return id2ann
    
def load_gtf_annotation(handle):
    """Parse sgd gff and return tuple of 3 elements:
    """
    id2ann = {}
    for line in handle:
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
        if feature!='gene': 
            continue
        id = comments.split('gene_id "')[-1].split('"')[0]
        note = comments.split('note "')[-1].split('"')[0]
        id2ann[id] = (note,) 
        
    return id2ann

def load_annotation( handle,dotStrip=0 ):
    """Return dictionary of annotation for every gene
    """
    id2ann = {}
    for l in handle:
        l=l.strip()
        if not l or l.startswith("#"):
            continue
        if not l.split("\t")[1:]:
            continue
        id = l.split("\t")[0]
        if dotStrip and "." in id:
            rDotIndex = id.rindex(".")
            id = id[:rDotIndex]
        id2ann[ id ] = l.split("\t")[1:]
    return id2ann

def load_cuffdiff(handle):
    """Load cuffdiff info
    """
    id2exp   = {}
    for line in handle: 
        if line.startswith("test_id"):
            continue
        #split line
        data = line.split('\t')
        #old cufflinks version
        locus = ""
        if len(data)==10: 
            #test_id  gene  locus  status  value_1  value_2  log(fold_change)  test_stat  p_value  significant
            transid,geneid,status,v1,v2,log_fc,test_stat,p,significant = data
        #ver >1.0.3
        else:             
            #test_id  gene_id  gene  locus  sample_1  sample_2  status  value_1  value_2  ln(fold_change)  test_stat  p_value  q_value  significant
            transid,geneid,gene,locus,s1,s2,status,v1,v2,log_fc,test_stat,p,q,significant = data
      
        v1,v2,log_fc,p=float(v1),float(v2),float(log_fc),float(p)
        id2exp[transid]=[gene,locus,v1,v2,log_fc,p] 
        
    return id2exp

def report(files, pfam, annotation, tab, pTh, verbose):
    """ """
    #load pfam annotation
    geneid2pfam, geneid2annotation, geneid2tab = {}, {}, {}
    if pfam:
        geneid2pfam = load_pfam_tblout(pfam)
        sys.stderr.write(" PFAMs for %s entries loaded.\n" % len(geneid2pfam))
    if tab:
        geneid2tab = load_annotation(open(tab))
        sys.stderr.write(" Tab annotation for %s entries loaded.\n" % len(geneid2tab))  
    if annotation:
        if annotation.endswith('.gff'):
            geneid2annotation = load_gff_annotation(open(annotation))
        elif annotation.endswith('.gtf'):
            geneid2annotation = load_gtf_annotation(open(annotation))
        else:
            geneid2annotation = load_annotation(open(annotation),True )
        sys.stderr.write( " Annotations for %s entries loaded.\n" % len(geneid2annotation) )
                        
    #load all cuffdiff files
    fn2data = {}
    fnames  = []
    for f in files:
        fn = f.name
        fnames.append(fn)
        fn2data[fn] = load_cuffdiff(f)

    #write output
    header = "#transcript id\tgene id"
    for fn in fnames:
        header += "\tcontrol\t%s\tlog2(FC)\tP-value" % fn
    header += "\tannotation\n"
    sys.stdout.write( header )
    #open outfiles for ids
    outfiles = [ open( "%s.%s.ids" % (fn,pTh),"w") for fn in fnames  ]
    #process all genes
    for transid in sorted(fn2data[fn]):
        lineData = []
        pFilter=False
        exprData = [ fn2data[fn][transid] for fn in fnames ] #; print exprData
        for exprTuple,out in zip(exprData,outfiles):
            geneid,locus,v1,v2,log_fc,p = exprTuple
            if not lineData:
                lineData = [ transid,geneid ]
            lineData += [ str(v1), str(v2), str(log_fc),str(p) ]
            #check p value
            if p<=pTh:
                pFilter = True
                out.write( geneid+"\n" )
                
        if pFilter:
            #add PFAM annotation
            annList=[]
            if transid in geneid2pfam:
                if type(geneid2pfam[transid]) is list:
                    annList.append(";".join(geneid2pfam[transid]))
                else:
                    for pfam,data in geneid2pfam[transid].iteritems():
                        annList.append( "%s [%s]" % (data[1],pfam) ) 
            lineData.append( "; ".join(annList))
            #add tab annotation
            annList=[]
            if geneid2tab and transid in geneid2tab:
                annList.append(";".join(geneid2tab[transid]))
            lineData.append( "; ".join(annList))
            #add Arabidopsis annotation
            if transid in geneid2annotation:
                for ann in geneid2annotation[transid]:
                    lineData.append( ann )
            #output info
            sys.stdout.write( "\t".join( lineData ) + "\n" )

    for out in outfiles:
        out.close()
    
def main():
    """
    """
    usage   = "usage: %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true")
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="files", nargs="+", type=file,
                        help="cuffdiff .diff file(s)   [%(default)s]")
    parser.add_argument("-p", dest="pTh", default=0.05, type=float,
                        help="P-value cut-off          [%(default)s]")
    parser.add_argument("-a", dest="annotation", default='',# type=file, 
                        help="annotation file          [%(default)s]" )
    parser.add_argument("-b", "--pfam", default='', #type=file,
                        help="pfam annotation file     [%(default)s]" )
    parser.add_argument("-t", "--tab", default='', #type=file,
                        help="pfam annotation file     [%(default)s]" )
  
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write( "Options: %s\n" % str(o) )

    report(o.files, o.pfam, o.annotation, o.tab, o.pTh, o.verbose)

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "#Time elapsed: %s\n" % dt )
