#!/bin/bash 

REF=$1  ## reference genome fasta 
STR=${REF%%_*}

cd fastqs 
KK=`ls | grep $STR | grep R1`
## this is done in case you have more than 1 strain 
## in our case, two strains (P125109 and D7795) were pro

for i in $KK
do
  TAG=${i%%.R1.fastq.gz}
  echo "BWA: Processing sample $TAG.." 
  R1=$i
  R2=$TAG.R2.fastq.gz
  bwa mem -t 16 $REF $R1 $R2 | samtools sort -@16 -O BAM - > $TAG.bam & 
done 
wait 

mv $STR*bam ../bams 

echo "BWA PAIRED-END ALIGNMENT: ALL DONE!"
