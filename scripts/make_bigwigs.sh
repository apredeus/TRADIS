#!/bin/bash 

## right utility allows you not to think about reference here

cd bams 

for i in *filt.bam
do
  echo "Indexing file $i.."
  samtools index $i & 
done 
wait

for i in *.filt.bam 
do
  TAG=${i%%.filt.bam}
  echo "deepTools: making bigWig file for sample $TAG.."
  bamCoverage -b $i -bs 1 -o $TAG.bw &> /dev/null &  
done 
wait 

mv *.bw ../browser

echo "MAKING BIGWIGS: ALL DONE!" 
