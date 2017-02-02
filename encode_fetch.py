#!/usr/bin/env python

#!/usr/bin/env python
desc="""
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Warsaw, 1/02/2017
"""

import os, sys, gzip
from datetime import datetime

import requests, json
 
# Force return from the server in JSON format
HEADERS = {'accept': 'application/json'}

def _get_response(url, fname="", jsondir=".json"):
    """Return server response"""
    if not os.path.isdir(jsondir):
        os.makedirs(jsondir)
    # load from file
    outfn = os.path.join(jsondir, "%s.json.gz"%fname)
    if fname and os.path.isfile(outfn):
        return json.load(gzip.open(outfn))
    # fetch 
    response = requests.get(url, headers=HEADERS)
    info = response.json()
    # dump                     
    if fname:
        with gzip.open(outfn, "w") as out:
            json.dump(info, out)
    return info

def get_experiments(assay_title="RNA-seq", www="https://www.encodeproject.org", end="&frame=object&format=json", verbose=1,
                    query="/search/?type=Experiment&assay_title=%s&limit=all&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens"):
    """Return experiments for given taxa""" #&award.rfa=ENCODE3
    url = "".join([www, query%assay_title, end])
    info = _get_response(url)
    experiments = [info['@graph'][i]['accession'] for i in range(len(info['@graph']))]
    if verbose:
        logger(" %s %s experiments found."%(assay_title, len(experiments)))
    return sorted(experiments)

def _get_donor_info(info):
    """Return dictionary of accession and donor info"""
    acc2donor = {}
    for f in info['files']:    
        acc = f["accession"]
        if not "replicate" in f:
            continue 
        # /human-donors/ENCDO845WKR/ -> ENCDO845WKR
        donor = f["replicate"]["library"]["biosample"]["donor"].split('/')[-2]
        name = f["replicate"]["library"]["biosample"]["biosample_term_name"]
        acc2donor[acc] = (donor, name)
    return acc2donor
    
def update_donors(donors, experiment, www="https://www.encodeproject.org/experiments/%s/?format=json"):
    """Return donors for given experiment"""
    # get json
    info = _get_response(www%experiment, experiment)
    # get donor info before as some exp lack donor info for not-fastq files
    acc2donor = _get_donor_info(info)
    # get requested files
    for f in info['files']:
        if f["status"] != "released" or f["file_type"] != "bam" or f["output_type"] != "alignments" \
           or f["assembly"] not in ["GRCh38", "hg19"]:
            continue
        acc, href = f["accession"], f["href"]
        # get parent / derived_from    
        pacc = None
        if "derived_from" in f:
            ids = set(f["derived_from"][i]["accession"] for i in range(len(f["derived_from"]))
                      if "accession" in f["derived_from"][i])
            if ids:
                pacc = ids.pop() #f["derived_from"][0]["accession"]
        # get donor info
        if acc in acc2donor:
            donor, name = acc2donor[acc]
        elif pacc and pacc in acc2donor:
            donor, name = acc2donor[pacc]
        else:
            continue
        # update donors
        if donor not in donors:
            donors[donor] = {}
        donors[donor][acc] = (name, href)
    return donors

def encode_fetch(out, assays, verbose=1):
    """Save experiments having RNase-seq and RNA-seq from the same donors"""
    if verbose:
        logger("Parsing experiments...")    
    donor2experiments = {}
    for assay_title in assays:
        donor2experiments[assay_title] = {}
        if not os.path.isdir(assay_title):
            os.makedirs(assay_title)
        for i, experiment in enumerate(get_experiments(assay_title), 1):
            sys.stderr.write("  %s %s donors\r"%(i, len(donor2experiments[assay_title])))
            donor2experiments[assay_title] = update_donors(donor2experiments[assay_title], experiment)
        logger("  %s donors for %s runs"%(len(donor2experiments[assay_title]), sum(len(x) for x in donor2experiments[assay_title].itervalues())))

    # print common
    intersection = set(donor2experiments[assays[0]]).intersection(donor2experiments[assays[1]])
    logger("Donors intersection: %s"%len(intersection), 1)

    www = "https://www.encodeproject.org"
    for donor in intersection:
        for a in assays:
            for acc, (name, href) in donor2experiments[a][donor].iteritems():
                name = name.replace(' ', '_').replace("'", '_')
                fname = "%s.%s.%s.bam"%(donor, name, acc)
                outfn = os.path.join(a, fname)
                if not os.path.isfile(outfn):
                    cmd = "wget -O %s -nc %s%s"%(outfn, www, href)
                    out.write(cmd+"\n")
        
def logger(info="", raw=0):
    sys.stderr.write(info+'\n')
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0a')   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")
    parser.add_argument("-a", "--assays", nargs="+", default=["RNA-seq", "DNase-seq"],
                        help="assay titles [%(default)s]")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")    
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    encode_fetch(o.output, o.assays, o.verbose)

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
