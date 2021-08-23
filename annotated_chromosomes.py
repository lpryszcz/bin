#!/usr/bin/env python2
desc="""Generate chromosome graphs for given genome.

http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec343

"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Barcelona, 13/05/2014
"""

import argparse, os, sys
from datetime import datetime
from reportlab.lib import colors
#from genome_annotation import get_contig2coverage,load_sgd_gff,parse_gtf

from Bio.SeqFeature import SeqFeature, FeatureLocation

from reportlab.lib.units import cm
from Bio import SeqIO
from Bio.Graphics import BasicChromosome

from genome_track import load_counts_beds, get_log2

def get_features(bed1, expcount1, bed2, expcount2, window):
    """Generate feature data:
    - mark heterozygous SNP-rich with orange
    - mark homozygous SNP-rich with darkgrey
    """
    def add_feature(d, gdata, limit):
        s,e,strand,name,color = d
        if e-s>=limit and color:
            gdata.append((s,e,strand,name,color))
        return gdata
                         
    gdata = []
    bed1.sort()
    bed2.sort()
    d=[0,0,0,0,0]
    limit = 1*window
    for (s,e,c1), (s2,e2,c2) in zip(bed1, bed2):
        #white for hetero
        if   get_log2(s,e,c2,expcount2,window)>0:
            if d[-1] == "":
                d[1] = e
            else:
                gdata = add_feature(d, gdata, limit)
                d = [s,e,0,'',''] #1
        #orange for hapB
        elif get_log2(s,e,c1,expcount1,window)>0:
            if d[-1] == "orange":
                d[1] = e
            else:
                gdata = add_feature(d, gdata, limit)
                d = [s,e,0,'','orange'] #1
        #grey for hapA
        else:
            if d[-1] == "grey":
                d[1] = e
            else:
                gdata = add_feature(d, gdata, limit)
                d = [s,e,0,'','grey'] #-1
        gdata = add_feature(d, gdata, limit) 
    return gdata

def get_features0(bed1, expcount1, bed2, expcount2, window):
    """Generate feature data:
    - mark heterozygous SNP-rich with orange
    - mark homozygous SNP-rich with darkgrey
    """        
    gdata = []
    bed1.sort()
    bed2.sort()
    for (s,e,c1), (s2,e2,c2) in zip(bed1, bed2):
        #white for hetero
        if   get_log2(s,e,c2,expcount2,window)>0:
            continue
        #orange for hapB
        elif get_log2(s,e,c1,expcount1,window)>0:
            gdata.append((s,e,1,'','orange'))
        #grey for hapA
        else:
            gdata.append((s,e,-1,'','grey'))
    return gdata
    

def annotated_chromosomes(fasta, output, spname, homosnps, heterosnps, scale, \
                          telomere_length, window, lenlimit, verbose, multi=10 ):
    """Generate chromosome plot"""
    #load bed files
    homocountsdict, expcounts1, homofns = load_counts_beds(homosnps, window, 0, verbose)
    hetecountsdict, expcounts2, hetefns = load_counts_beds(heterosnps, window, 0, verbose)
    expcount1 = expcounts1[0]
    expcount2 = expcounts2[0]
    #get chromosome names and lengths
    chr2length = {r.id: len(r) for r in SeqIO.parse(fasta, 'fasta')}
    #total genome length
    max_len = max(chr2length.values())
    if verbose:
        sys.stderr.write("%s chromosomes. The largest chromosome is %s bp\n"%(len(chr2length), max_len))
    #init diagram
    chr_diagram = BasicChromosome.Organism()
    multisize = 5
    chr_diagram.page_size = (multi*29.7*cm*multisize, multi*21*cm*multisize) #A4 landscape
    chr_diagram.output_format=output.split('.')[-1]
    chr_diagram.title_size=20*multi
    #add chromosomes
    for i, (name, length) in enumerate(sorted(chr2length.items(), key=lambda x: x[1], reverse=True)):
        '''features = [f for f in record.features if f.type=="tRNA"]
        #Record an Artemis style integer color in the feature's qualifiers,
        #1 = Black, 2 = Red, 3 = Green, 4 = blue, 5 =cyan, 6 = purple 
        for f in features: f.qualifiers["color"] = [index+2]'''
        if length<lenlimit*1e3:
            continue
        print i, name, length
        cur_chromosome = BasicChromosome.Chromosome(name.split()[0].split('|')[0])
        #Set the scale to the MAXIMUM length plus the two telomeres in bp,
        #want the same scale used on all five chromosomes so they can be
        #compared to each other
        cur_chromosome.scale_num = max_len + 2 * telomere_length
        cur_chromosome.title_size = 12*multi
            
        #Add an opening telomere
        start = BasicChromosome.TelomereSegment()
        start.scale = telomere_length
        cur_chromosome.add(start)
        
        #get counts
        bed1 = bed2 = ([],[]), ([],[])
        if homocountsdict:
            bed1 = homocountsdict[name][0]
        if hetecountsdict:
            bed2 = hetecountsdict[name][0]
        features = get_features(bed1, expcount1, bed2, expcount2, window)
        #add scale marker
        if not i:
            for i in xrange(0, length, int(scale/2)):
                features.append((i, i+1, 0, "%.2f Mb"%(i/scale,), 'black'))
            
        #Add a body - again using bp as the scale length here.
        body = BasicChromosome.AnnotatedChromosomeSegment(length, features)
        body.scale = length
        body.label_size = 6*multi
        cur_chromosome.add(body)
        
        #Add a closing telomere
        end = BasicChromosome.TelomereSegment(inverted=True)
        end.scale = telomere_length
        cur_chromosome.add(end)

        #This chromosome is done
        chr_diagram.add(cur_chromosome)

    #draw
    chr_diagram.draw(output, spname)

def main():
    #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog)

    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')
    parser.add_argument("-i", "--input", required=True,
                        help="genome fasta file" )
    parser.add_argument("-o", "--output", required=True, 
                        help="output file name" )
    parser.add_argument("-n", "--name", default="", 
                        help="Species name [%(default)s]" )
    parser.add_argument("-s", dest="homosnp",  nargs='+', type=file,
                        help="Homozygous SNP BED files [%(default)s]" )
    parser.add_argument("-t", dest="hetesnp",  nargs='+', type=file,
                        help="Heterozygous SNP BEDs files [%(default)s]" )    
    parser.add_argument("-w", dest="window",   default=1000, type=int, 
                        help="window size         [%(default)s bp]" )
    parser.add_argument("-l", dest="lenlimit", default=100, type=int, 
                        help="only chr above l kb [%(default)s kb]" )
    parser.add_argument("--telomere_length", default=3*1e4, 
                        help="default length of telomeres [%(default)s]" )
    parser.add_argument("--chromosome_scale", default=1e6, 
                        help="chromosome scale [%(default)s]" )

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    annotated_chromosomes(o.input, o.output, o.name, o.homosnp, o.hetesnp, o.chromosome_scale, \
                          o.telomere_length, o.window, o.lenlimit, o.verbose)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
