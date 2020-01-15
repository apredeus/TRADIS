#!/bin/bash

## conda activate tradis - see setup instruction on https://github.com/apredeus/TRADIS

R1="ls raw/* | tr '\n' '\t' | cut -f 1"
R2="ls raw/* | tr '\n' '\t' | cut -f 2"
BC="barcodes.fa" 

cutadapt -O 34 -g file:$BC $R1 $R2 -o fastqs/{name}.R1.fastq.gz -p fastqs/{name}.R2.fastq.gz --discard-untrimmed &> cutadapt.log 
