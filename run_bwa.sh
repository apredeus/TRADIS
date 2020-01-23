#!/bin/bash 

STR=$1  ## reference genome should be named STR.fa 
PREFIX=${STR%%_*}

cd fastqs 
KK=`ls | grep $PREFIX | grep R1`
## this is done in case you have more than 1 strain 
## in our case, two strains (P125109 and D7795) were processed in parallel

for i in $KK
do
  TAG=${i%%.R1.fastq.gz}
  echo "BWA: Processing sample $TAG.." 
  R1=$i
  R2=$TAG.R2.fastq.gz
  bwa mem -t 16 $REF $R1 $R2 | samtools sort -@16 -O BAM - > $TAG.bam &> $TAG.bwa.log & 
done 
wait 

mv $PREFIX*bam ../bams 
mv $PREFIX*log ../logs 

echo "BWA PAIRED-END ALIGNMENT: ALL DONE!"
