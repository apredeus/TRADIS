#!/bin/bash 

TAG=$1
REF=$2

R1=fastqs/$TAG.R1.fastq.gz
R2=fastqs/$TAG.R2.fastq.gz

bwa mem -t 16 $REF $R1 $R2 | samtools sort -O BAM - > $TAG.bam
samtools index $TAG.bam
picard MarkDuplicates I=$TAG.bam O=$TAG.rmdup.bam M=$TAG.metrics.out REMOVE_DUPLICATES=true &> $TAG.rmdup.log 
