#!/bin/bash 

REF=$1
STR=${REF%%_*}

cd fastqs 
KK=`ls | grep $STR | grep R1`

echo "Processing samples: " 
echo $KK 

cd ..
 
for i in $KK
do
  TAG=${i%%.R1.fastq.gz}
  ./run_bwa.sh $TAG $REF &> $TAG.bwa.log &  
done 

wait 

mv $STR*bam $STR*bai ../bams 
mv $STR*log $STR*metrics.out ../logs

echo "ALIGNMENT AND DUPLICATE REMOVAL: ALL DONE!"
