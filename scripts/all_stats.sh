#!/bin/bash 

## calculate all alignment stats etc

cd bams 
KK=`ls *.filt.bam | sed "s/.filt.bam//g"`
cd ..

echo -e "Sample\tN_reads\tN_aligned\tPct_aligned\tN_dedup\tPct_dup\tN_filt\tN_gene\tPct_gene\tN_uniq\tN_ugene\tPct_ugene"

for i in $KK
do
  ./calculate_stats.sh $i & 
done 

wait 
