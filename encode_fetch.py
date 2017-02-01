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

def _get_response(url, fname=""):
    """Return server response"""
    # load from file
    outfn = "%s.json.gz"%fname
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
                    query="/search/?type=Experiment&assay_title=%s&limit=all&award.rfa=ENCODE3&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens"):
    """Return experiments for given taxa"""
    url = "".join([www, query%assay_title, end])
    info = _get_response(url)
    experiments = [info['@graph'][i]['accession'] for i in range(len(info['@graph']))]
    if verbose:
        logger(" %s %s experiments found."%(assay_title, len(experiments)))
    return sorted(experiments)

def update_donors(donors, experiment, www="https://www.encodeproject.org/experiments/%s/?format=json"):
    """Return donors for given experiment"""
    info = _get_response(www%experiment, experiment)
    for f in info['files']:
        if f["status"] != "released" or f["file_type"] != "bam" or f["output_type"] != "alignments":
            continue
        acc, href = f["accession"], f["href"]
        if not "replicate" in f:
            #logger("[ERROR] Cannot process %s"%experiment) #print experiment, f
            continue #sys.exit(1)
        # /human-donors/ENCDO845WKR/ -> ENCDO845WKR
        donor = f["replicate"]["library"]["biosample"]["donor"].split('/')[-2]
        name = f["replicate"]["library"]["biosample"]["biosample_term_name"]
        # update donors
        if donor not in donors:
            donors[donor] = {}
        donors[donor][acc] = (name, href)
    return donors

def encode_fetch(assays, verbose=1):
    """Save experiments having RNase-seq and RNA-seq from the same donors"""
    if verbose:
        logger("Parsing experiments...")    
    donor2experiments = {}
    for assay_title in assays:
        donor2experiments[assay_title] = {}
        for experiment in get_experiments(assay_title):
            donor2experiments[assay_title] = update_donors(donor2experiments[assay_title], experiment)

    # print summary
    info = "Donors:" 
    for a in assays:
        info += "\n %s for %s"%(len(donor2experiments[a]), a)
    # common
    intersection = set(donor2experiments[assays[0]]).intersection(donor2experiments[assays[1]])
    info += "\n%s for both"%len(intersection)
    logger(info, 1)
        
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
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    encode_fetch(o.assays, o.verbose)

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
