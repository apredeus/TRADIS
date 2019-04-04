#!/bin/bash 

FA=$1
samtools faidx $FA
cut -f 1,2 $FA.fai > ${FA%%.fa}.chrom.sizes

STR=${FA%%_*}
CHROM=`readlink -f ${FA%%.fa}.chrom.sizes`

cd bams 
KK=`ls *.filt.bam | grep $STR`

for i in $KK
do
  echo "Processing sample $i"
  bedtools genomecov -ibam $i -bg > ${i%%.filt.bam}.bedGraph & 
done 
wait

for i in *.bedGraph
do
  bedGraphToBigWig $i $CHROM ${i%%.bedGraph}.bw &
done 
wait 

rm *.bedGraph
echo "MAKING BIGWIGS: ALL DONE!" 
