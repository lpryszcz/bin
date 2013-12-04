#!/usr/bin/env python
###
#
###

import os,sys
from datetime import datetime
from Bio import SeqIO

def main():
  
  for i, r in enumerate(SeqIO.parse(sys.stdin, 'fasta'), 1):
  	contig = "contig%5i" % i
  	contig = contig.replace(' ', '0')
  	sys.stdout.write(">%s\n%s\n" % (contig, "\n".join(str(r.seq[j:j+60]) for j in range(0, len(r), 60))))
  	'''#r.description
    header = r.id
    
    contig="contig%s" % i
    if not header.startswith( ('supercontig','contig') ): 
      while contig in contigs: 
        i+=1
        contig="contig%s" % i
      header=contig #"%s lenght:%s" % ( contig,len(r.seq) )
      contigs.append(contig)
    else: contigs.append(r.id)
    
    j=0
    fasta=">%s\n%s\n" % ( header,r.seq[j:j+window])
    while j<=len(r.seq):
      j+=window
      fasta+='%s\n' % r.seq[j:j+window]
      
    out.write( fasta )
    
    n=str(r.seq).count('N')
    sys.stderr.write( "%s %s Ns (%2.2f%s)\n" % (header,n,n*1.0/len(r.seq),'%') )
  out.close()'''

if __name__=='__main__': 
  t0=datetime.now()
  main()
  dt=datetime.now()-t0
  sys.stderr.write( "Time elapsed: %s\n" % dt )
