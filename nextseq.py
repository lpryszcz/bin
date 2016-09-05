#!/usr/bin/env python
desc="""Converts binary Illumina files (.bcl) to FastQ. Designed to work as cron service.
It proceeds only if given run has been sequenced (.bcl present),
but not yet successfully converted (no .fastq files). 

Dependencies:
- bcl2fastq 
"""
epilog="""Author: l.p.pryszcz@gmail.com
Warsaw, 11/01/2016
"""

import os, sys, subprocess, re
from datetime import datetime
from glob import glob
import xml.etree.ElementTree as ET

def _isCompleted(fpath):
    """Return True if run completed"""
    xml = ET.parse(open(fpath))
    e = xml.find('CompletionStatus')
    if e.text == "CompletedAsPlanned":
        return True
    else:
        sys.stderr.write("[WARNING][isCompleted] %s says %s\n"%(fpath, e.text))

def _isProcessed(outdir):
    """Return True if current run is successfully processed.
    If there are some files missing, remove output directory
    and rerun conversion.
    """
    reportFn = os.path.join(outdir, "Reports/html/index.html")
    # not output directory
    if not os.path.isdir(outdir):
        return
    # missing report - uncompleted run
    elif not os.path.isfile(reportFn):
        # remove outdir
        sys.stderr.write("[WARNING] Incomplete conversion detected: %s. Remove it!\n"%outdir) #Cleaning up & reruning...
        #os.system("rm -r %s"%outdir)
        #return
    return True
        
def get_experiment_id(fpath):
    """Return True if run not processed"""
    xml = ET.parse(open(fpath))
    # get experiment name and flowcell id
    e = xml.find('ExperimentName') # LibraryID # 
    expName = e.text
    e = xml.find('FlowCellSerial')
    flowcellID = e.text
    #
    expID = "%s.%s" % (expName, flowcellID)
    # return safe filename # http://stackoverflow.com/a/295147/632242
    expID = re.sub('[^-a-zA-Z0-9_.() ]+', '', expID).replace(" ","_")
    return expID

def get_new_runs(fqdir, raw):
    """Return finished runs that are not yet processed."""
    for fpath in sorted(glob(os.path.join(raw, "*", "RunCompletionStatus.xml"))):
        rundir = os.path.dirname(fpath)
        runID  = os.path.basename(rundir)
        outdir = os.path.join(fqdir, runID)
        sampleSheetFpath = os.path.join(rundir, "SampleSheet.csv")
        # report only if run completed and given runID is not already processed AND SampleSheet.csv present
        if _isCompleted(fpath) and not _isProcessed(outdir) and os.path.isfile(sampleSheetFpath):
            fpath2 = os.path.join(rundir, "RunParameters.xml")
            expID = get_experiment_id(fpath2)
            yield rundir, outdir, runID, expID
        
def process_runs(fqdir, raw, threads, bindir, verbose, bufsize=-1):
    """Convert binary NextSeq files (.bcl) to FastQ."""
    # process unconverted runs
    for i, (curdir, outdir, runID, expID) in enumerate(get_new_runs(fqdir, raw), 1):
        if verbose:
            sys.stderr.write("[%s] %s %s\n"%(datetime.ctime(datetime.now()), runID, expID))
        # bcl2fastq
        #prepare outdir and paths
        #outdir = os.path.join(fqdir, runID)
        interopdir = os.path.join(outdir, "InterOp")
        os.makedirs(outdir)
        # define args
        threads = str(threads)
        args = ["%sbcl2fastq"%bindir, "--runfolder-dir", curdir, "--output-dir", outdir,
                "--interop-dir", interopdir, "--no-lane-splitting",
                "-r", threads, "-d", threads, "-p", threads, "-w", threads]
        #execute
        with open(os.path.join(outdir, "bcl2fastq.log"), "w") as logFile:
            proc = subprocess.Popen(args, bufsize=bufsize, stdout=logFile, stderr=logFile)
            proc.wait()
        # link experiment to run id
        cmd = "ln -s %s %s"%(runID, os.path.join(fqdir, expID))
        os.system(cmd)

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-d", "--dir", default="/mnt/illumina_fastq/lpryszcz/nextseq",
                        help="working (.fastq) directory [%(default)s]")
    parser.add_argument("-r", "--raw", default="/mnt/illumina_raw/NextSeq_data",
                        help="raw (.bcl) directory       [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int, 
                        help="threads                    [%(default)s]")
    parser.add_argument("--bindir", default="",
                        help="binary directory prefix    [%(default)s]")
   
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    process_runs(o.dir, o.raw, o.threads, o.bindir, o.verbose)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("%s\n"%str(e)) #"I/O error({0}): {1}\n".format(str(e), e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)

