#!/bin/bash 

## DE_GFF - GFF file where every feature is a "gene" and ID=locus tag
## ESS_GFF - same, but with last 10% of the gene removed (equivalent of -trim3 0.1) 

TAG=$1              ## e.g. P125109_v3_2018
STR=${TAG%%_*}      ## e.g. P125109 
DIR=`pwd` 
DE_GFF=$DIR/$TAG.tradis.gff
ESS_GFF=$DIR/$TAG.ess_90.gff

cd bams 
KK=`ls *deseq.bam | grep $STR`
PP=`ls *ess.bam   | grep $STR`

for i in $KK
do
  echo "featureCounts: processing DESeq2 sample $i.." 
  TAG=${i%%.bam}
  ## no strand-specificity, multimappers/multifeatures are split proportionally 
  featureCounts -M -O --fraction -t gene -g ID -s 0 -a $DE_GFF -o $TAG.fc.tsv $i &> $TAG.fc.log &  
done 
wait

for j in $PP
do
  echo "featureCounts: processing essentiality sample $j.." 
  TAG=${j%%.bam}
  ## no strand-specificity, multimappers are simply counted like normal reads, and all overlapping features are counted twice
  featureCounts -M -O -t gene -g ID -s 0 -a $ESS_GFF -o $TAG.fc.tsv $j &> $TAG.fc.log &  
done 
wait

rm *summary 
mv *.fc.tsv ../counts
mv *.fc.log ../logs 

echo "FEATURE COUNTS: ALL DONE!"
