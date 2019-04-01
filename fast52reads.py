#!/usr/bin/env python
desc="""

TDB:
- 

CHANGELOG:
v1.1
- 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Warsaw, 18/06/2018
"""

import os, sys
from datetime import datetime
import numpy as np
#from multiprocessing import Pool
#from sklearn.preprocessing import normalize
#from scipy import stats, signal
import matplotlib.pyplot as plt
import h5py

true, false, null = True, False, False

def get_metadata(handle):
    """Return metadata"""
    exec("meta="+handle['Meta']['User']['analysis_conf'][0])
    # meta['experiment']: ['sequencing', 'GUI', 'flowcell', 'mux_scan', 'kit']
    # meta['channel_states']['3'].keys() ['style', 'group', 'description', 'name', 'logic']
    # meta['histograms'].keys() ['read_length', 'read_median', 'event_length', 'read_event_count_weighted', 'event_mean']
    # get channel names
    meta["channels"] = list(handle['Raw'].keys())
    # sample rate +  offset + range (range_in_pA)
    meta["sample_rate"] = handle['IntermediateData'][meta["channels"][0]]['Meta'].attrs["sample_rate"]
    return meta

def get_active_pores(handle, min_pA=10, nth=100, logger=0):
    """Return the dictionary of active pores with name of the pore as key and mean and std as values."""
    #active = {ch: (np.mean(handle['Raw'][ch]['Signal'][::nth]), np.std(handle['Raw'][ch]['Signal'][::nth]))
    #          for ch in handle['Raw'].keys() if np.mean(handle['Raw'][ch]['Signal'][::nth])>min_pA}
    active = {}
    for ch in handle['Raw'].keys():
        sig = np.array(handle['Raw'][ch]['Signal'])
        if sig.mean()>min_pA:
            active[ch] = (sig.mean(), sig.std())
    return active

def get_reads(out, ch, sig, avg, stdev, minlen=100, logger=0):
    """Return reads from squiggle"""
    i = 0
    # get beginning of the read
    start = np.argmax(sig>avg)
    while start:
        th = avg - stdev
        # get start of the read
        read_start = np.argmax(sig[start:]<th)+start
        # break if no new read
        if read_start==start:
            break
        # get end of the read and its length
        read_end = np.argmax(sig[read_start:]>avg)+read_start
        rlen = read_end-read_start
        # set new start
        start = read_end
        # skip if too short read
        if rlen<minlen:
            continue
        i += 1
        info = "%s\t%s\t%s\t%s\t%s"%(ch, i, read_start, read_end, rlen)
        logger("%s\n"%info)
        #offset = 10; plt.plot(sig); plt.xlim(read_start-offset,read_end+offset); plt.ylim(0, 400); plt.show()
        # save read
        if out:
            out.write("%s\t%s\n"%(info, ";".join(map(str, sig[read_start:read_end]))))
    
def fast52reads(fast5, out, minlen, logger=0):
    """Process raw fast5 and report reads"""
    # open hdf file
    handle = h5py.File(fast5)

    # load metadata
    meta = get_metadata(handle)
    fc, kit = meta['experiment']['flowcell'], meta['experiment']['kit']
    fc_info = "Flow-cell: %s / Kit: %s\n"%(fc, kit)
    logger(fc_info)
    
    # get active pores
    active = {u'Channel_4': (234.41526257623752, 122.36777796984582), u'Channel_66': (275.32927299536664, 104.2860897951875)}#, u'Channel_140': (1024.8175093671694, 1158.4218179410236), u'Channel_245': (257.05529806510333, 99.62405722286907), u'Channel_158': (309.387211671434, 333.3532719811441), u'Channel_339': (304.4667964132901, 372.11462437007214), u'Channel_316': (415.7896542273415, 505.4421564545962), u'Channel_287': (298.0238700297858, 119.58053404885736), u'Channel_164': (469.3899645702331, 744.9271750383655), u'Channel_280': (202.92503826651222, 137.06719825641107), u'Channel_189': (221.0037620118434, 118.0202397377132), u'Channel_204': (362.206857831781, 379.241978230955), u'Channel_172': (95.64020502872205, 400.6039203998252), u'Channel_199': (288.90700314110444, 559.3561514336901), u'Channel_380': (229.4527423348778, 117.28768838963842), u'Channel_450': (288.54346241016975, 119.91383950600196), u'Channel_117': (160.67694014171906, 189.7024130745043), u'Channel_116': (628.0186693300553, 852.5172231443868)}
    active = get_active_pores(handle)
    pore_info = "%s active pores: %s\n"%(len(active), "\n".join(" %s: %.2f +- %.f"%(k, v[0], v[1]) for k, v in active.items()))
    logger(pore_info)

    if out:
        out.write("##%s\n##%s"%(fast5, fc_info))
        out.write("##%s"%pore_info.replace('\n', '\n##'))
        out.write("\n#channel\tread\tstart\tend\tlength\tpA readouts separated by ;\n")
    
    # get reads
    for channel in active:
        avg, stdev = active[channel]
        sig = np.array(handle['Raw'][channel]['Signal'])
        get_reads(out, channel, sig, avg, stdev, logger=logger)
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-i", "--fast5", nargs="+", help="BAM file")
    parser.add_argument("-o", "--output", default="", type=argparse.FileType('w'), help="output stream [stdout]")
    parser.add_argument("-l", "--minlen", default=100, type=int, help="min read length [%(default)s]")
    parser.add_argument("--log", default=sys.stderr, type=argparse.FileType('w'), help="redirect verbose to [stderr]")
    
    o = parser.parse_args()
    if o.verbose:
        logger = lambda x: o.log.write(x)
        logger("Options: %s\n"%str(o))
    else:
        logger = lambda x: not x

    # calculate coverage from bam for given intervals #get_reads
    for fast5 in o.fast5:
        logger("Processing %s ...\n"%fast5)
        fast52reads(fast5, o.output, o.minlen, logger)
 
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

