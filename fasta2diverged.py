#!/usr/bin/env python
desc="""Return sequences diverged by given percentage compares to input.

LOH sizes are drawn from negative binomial distribution, given r (--nb_r) and p (--nb_p).
Default r/p are learned from C. orthopsilosis MCO456 genomic data (PMID: )
using negative binomial fit function from: http://bit.ly/1fL8dC5 .

Is negative binomial really the best for that?
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow, 4/03/2014
"""

import os, random, sys
from datetime import datetime
from Bio import SeqIO, Seq
import numpy as np

aminos = 'ACDEFGHIKLMNPQRSTVWY'
nucleotides = 'ACGT'

import scipy.special as special
import scipy.optimize as optimize
import numpy as np
import mpmath

class negBin(object):
    """
    http://www.nehalemlabs.net/prototype/blog/2013/11/11/negative-binomial-with-continuous-parameters-in-python/
    """
    def __init__(self, p = 0.1, r = 10):
        nbin_mpmath = lambda k, p, r: mpmath.gamma(k + r)/(mpmath.gamma(k+1)*mpmath.gamma(r))*np.power(1-p, r)*np.power(p, k)
        self.nbin = np.frompyfunc(nbin_mpmath, 3, 1)
        self.p = p
        self.r = r

    def mleFun(self, par, data, sm):
        '''
        Objective function for MLE estimate according to
        https://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation

        Keywords:
        data -- the points to be fit
        sm -- \sum data / len(data)
        '''
        p = par[0]
        r = par[1]
        n = len(data)
        f0 = sm/(r+sm)-p
        f1 = np.sum(special.psi(data+r)) - n*special.psi(r) + n*np.log(r/(r+sm))
        return np.array([f0, f1])

    def fit(self, data, p = None, r = None):
        if p is None or r is None:
            av = np.average(data)
            va = np.var(data)
            r = (av*av)/(va-av)
            p = (va-av)/(va)
        sm = np.sum(data)/len(data)
        x = optimize.fsolve(self.mleFun, np.array([p, r]), args=(data, sm))
        self.p = x[0]
        self.r = x[1]

    def pdf(self, k):
        return self.nbin(k, self.p, self.r).astype('float64')

def seq2diverged(seq, divergence, loh, minsize, lohSizes, aminos, verbose):
    """Return diverged sequence"""
    seqlist = list(seq)
    #number of position to change
    k = int(round(divergence*len(seqlist)))
    positions = random.sample(xrange(len(seqlist)), k)
    #inclue LOHs
    if loh:
        lohs = []
        while len(positions) > (1-loh) * len(seq) * divergence:
            #get LOH random start
            s = random.randint(0, len(seq)-1)
            #and random LOH length from negative binomial distribution
            lsize = 0
            while lsize<minsize:
                lsize = round(random.sample(lohSizes, 1)[0])
            e = s + lsize
            lohs.append(lsize)
            #filter snp posiitons
            positions = filter(lambda p: p<s or p>e, positions)
        if verbose:
            sys.stderr.write(" %s LOHs (%s - %s) incorporated\n" % \
                             (len(lohs), min(lohs), max(lohs)))
    #change positions
    for i in positions:
        nb = seqlist[i]
        while nb == seqlist[i]:
            nb = random.choice(aminos)
        seqlist[i] = nb
    return "".join(seqlist)

def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-i", "--input",     default=sys.stdin, type=file, 
                        help="fasta file(s)   [stdin]")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-d", "--divergence",    default=0.01, type=float, 
                        help="divergence      [%(default)s]")
    parser.add_argument("--dna",             default=False, action='store_true',
                        help="DNA alphabet    [amino acids]")
    parser.add_argument("--loh",             default=0.00, type=float, 
                        help="level of LOH    [%(default)s]")
    parser.add_argument("--minsize",         default=100, type=int, 
                        help="min LOH size    [%(default)s]")
    parser.add_argument("--learn",    default=False, type=file, 
                        help="BED file to learn LOH sizes [%(default)s]")
    parser.add_argument("--nb_r",            default=0.36676450795532967, type=float, 
                        help="LOH negative binomial distribution r [%(default)s]")
    parser.add_argument("--nb_p",            default=0.9998294827977795, type=float, 
                        help="LOH negative binomial distribution p [%(default)s]")
    parser.add_argument("--power_a",         default=0.1, type=float, 
                        help="LOH power distribution a [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
    
    '''r = o.nb_r
    p = o.nb_p
    if o.learn:
        if o.verbose:
            sys.stderr.write("Fitting negative binomial to data...\n")
        #load LOH sizes from BED
        data = np.array([float(l.split()[2])-float(l.split()[1]) for l in o.learn])
        #fit data
        n = negBin()
        n.fit(data)
        #store
        r = n.r
        p = n.p
        if o.verbose:
            sys.stderr.write(" r=%f\n p=%f\n"%(r, p))

    #get LOH sizes distibution
    #lohSizes = np.random.negative_binomial(p, r, 100000)*100'''
    lohSizes = np.random.power(o.power_a, 1e5)*1e6
    if o.verbose:
        sys.stderr.write("Processing chromosomes...\n")
    for i, r in enumerate(SeqIO.parse(o.input, 'fasta'), 1):
        if o.verbose:
            sys.stderr.write('%s %s %sbp %s\n'%(i, r.id, len(r), " "*20))
        if o.dna:
            alphabet = nucleotides
        else:
            alphabet = aminos
        seq = seq2diverged(r.seq, o.divergence, o.loh, o.minsize, lohSizes, alphabet, \
                           o.verbose)
        r.seq = Seq.Seq(seq)
        o.output.write(r.format('fasta'))

if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
