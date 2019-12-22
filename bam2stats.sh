#!/bin/bash
# Report stats (mapped reads and identity to reference) from samtools stats 
# for bam file(s) ignoring secondary, suplementary and qc failed alignments
#
# USAGE: bam2stats.py bam1 bam2 ... bamN

for f in "$@"; do if [ -s $f ]; then echo $f `samtools stats -F3840 $f | awk '{if ($2=="sequences:"){nseq=$3} else if($2=="reads" && $3=="mapped:"){printf("mapped: %s (%.1f%)\n",$4,100*$4/nseq)} else if($2=="error" && $3=="rate:"){printf("identity: %.2f%\n",100-$4*1e2)}}'`
fi; done


