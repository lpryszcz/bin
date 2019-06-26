#!/bin/bash
# Report stats from samtools stats for bam file ignoring secondary, suplementary and qc failed alignments

samtools stats -F3840 $1 | awk '{if ($2=="sequences:"){nseq=$3} else if($2=="reads" && $3=="mapped:"){printf("mapped: %s (%.1f%)\n",$4,100*$4/nseq)} else if($2=="error" && $3=="rate:"){printf("identity: %.2f%\n",100-$4*1e2)}}'


