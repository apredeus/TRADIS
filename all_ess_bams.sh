#!/bin/bash 
#
# to give essentiality analysis more power, we need to merge BAM files and then find unique insertions

cd bams 
cp ../*.pl .

KK=`for i in *rmdup.bam ; do echo ${i%%_?.rmdup.bam}; done | sort | uniq`

for i in $KK
do
  echo "Processing sample $i (merging files BAM files).."
  samtools merge $i.merged.bam ${i}_?.bam & 
done 
wait 

for i in $KK
do
  echo "Processing sample $i (making uniq 1nt files).."
  samtools view -h -f 64 -F 4 $i.merged.bam | ./cigar_filter.pl | ./make_1nt_uniq.pl | samtools view -b - > $i.ess.bam & 
done
wait

echo "MERGED BAM PROCESSING: ALL DONE!" 
