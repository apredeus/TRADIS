#!/bin/bash 

TAG=$1
TAG2=${TAG%%_?} 

## all in bams folder 

N1=`samtools flagstat bams/$TAG.bam | grep read1 | awk '{print $1}'`  ## all read pairs
N2=`samtools view -f 64 -F 4 bams/$TAG.bam | wc -l`                   ## aligned
N3=`samtools view -f 64 bams/$TAG.rmdup.bam | wc -l`                  ## dedup 
N4=`samtools view bams/$TAG.filt.bam | wc -l`                         ## filtered/DEseq
N5=`samtools view bams/$TAG2.ess.bam | wc -l`                         ## unique insertions 

C1=`awk '{if (NR>2) sum+=$7} END {print sum}' counts/$TAG.deseq.fc.tsv`  
C2=`awk '{if (NR>2) sum+=$7} END {print sum}' counts/$TAG2.ess.fc.tsv`  

P1=`echo $N1 | awk -v v=$N2 '{printf "%.2f",v*100/$1}'`      ## % aligned 
P2=`echo $N2 | awk -v v=$N3 '{printf "%.2f",($1-v)*100/$1}'` ## % PCR duplicates 
P3=`echo $N4 | awk -v v=$C1 '{printf "%.2f",v*100/$1}'`      ## % filtered reads assigned to genes 
P4=`echo $N5 | awk -v v=$C2 '{printf "%.2f",v*100/$1}'`      ## % unique insertions assigned to genes 

echo -e "$TAG\t$N1\t$N2\t$P1\t$N3\t$P2\t$N4\t$C1\t$P3\t$N5\t$C2\t$P4" 
