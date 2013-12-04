#!/usr/bin/env python
desc="""Report histogram for given insert size data

"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 18/10/2012
"""

import argparse, os, sys
from datetime import datetime

def plot( isizes,outfn ):
    """
    """
    import matplotlib.pyplot as plt

    # the histogram of the data
    n, bins, patches = plt.hist(isizes, 50, normed=0, facecolor='g', alpha=0.75)

    plt.xlabel('Insert size')
    plt.ylabel('Occurencies')
    plt.title('Histogram of insert size [%s]' % os.path.basename(outfn).split(".")[0] )
    #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    #plt.axis([40, 160, 0, 0.03])
    plt.grid(True)

    if not outfn:
        plt.show()
    else:
        plt.savefig(outfn)

def sam2hist( handle,outfn,colnumb,verbose ):
    """
    """
    #
    isizes = []
    for l in handle:
        try:
            i = int(l.split('\t')[colnumb])
            isizes.append(i)
        except:
            if verbose:
                sys.stderr.write( "Warning: Cannot read isize in column %s for line: %s\n" % (colnumb,str(l.split('\t'))) )

    #plot
    plot( isizes,outfn )

def main():
    usage   = "usage: samtools view -f35 BAM [region] | %(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )

    parser.add_argument("-v", dest="verbose",  default=False, action="store_true")
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", dest="input",    default=sys.stdin, type=file, #argparse.FileType('r'),  
                        help="input sam file   [%(default)s]" )
    parser.add_argument("-c", dest="column",   default=8, type=int, 
                        help="column number 0-based [%(default)s]" )
    parser.add_argument("-o", dest="outfn",   default="", type=str, 
                        help="output fname     [%(default)s]" )    
    
    o = parser.parse_args()
    
    if o.verbose: 
        sys.stderr.write( "Options: %s\n" % str(o) )

    #create outdir
    outdir = os.path.dirname( o.outfn )
    if outdir and not os.path.isdir( outdir ):
        os.makedirs( outdir )
        
    sam2hist( o.input,o.outfn,o.column,o.verbose )


if __name__=='__main__': 
    t0  = datetime.now()
    main()
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )