#!/bin/bash 

cd bams 
cp ../*.pl .

for i in *.rmdup.bam 
do
  echo "Processing sample $i (1nt BAM files for featureCounts/DESeq2 processing).."
  TAG=${i%%.rmdup.bam} 
  ## here's what's happening here: 
  ## - samtools only keeps aligned R1 reads; 
  ## - cigar_filter.pl removes all reads that don't have salmonella seq right after the tag; 
  ## - make_1nt.pl truncates the BAM file to 1nt reads 

  samtools view -h -f 64 -F 4 $i | ./cigar_filter.pl | ./make_1nt.pl | samtools view -b - > $TAG.deseq.bam & 
done 
wait 

for i in *.rmdup.bam 
do
  echo "Processing sample $i (R1 BAM files for JBrowse visualization).."
  TAG=${i%%.rmdup.bam} 
  ## here's what's happening here: 
  ## - samtools only keeps aligned R1 reads; 
  ## - cigar_filter.pl removes all reads that don't have salmonella seq right after the tag; 

  samtools view -h -f 64 -F 4 $i | ./cigar_filter.pl | samtools view -b - > $TAG.jbrowse.bam & 
done 
wait 

echo "JBROWSE/DESEQ2 BAM PROCESSING: ALL DONE!" 
