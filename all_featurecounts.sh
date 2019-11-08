#!/bin/bash 

## this should be a genes file with ID=locus tag
## this way you can include both CDS and ncRNA 
GFF=$1
STR=${GFF%%_*} 
GFF=`readlink -f $GFF`

cd bams 
KK=`ls *deseq.bam *ess.bam | grep $STR`

for i in $KK
do
  echo "Processing sample $i" 
  TAG=${i%%.bam}
  featureCounts -M --fraction -t gene -g ID -s 0 -a $GFF -o $TAG.fc.tsv $i &> $TAG.fc.log &  
done 
wait

rm *summary 
mv *.fc.tsv *.fc.log ../counts 

echo "FEATURE COUNTS: ALL DONE!"
