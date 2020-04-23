#!/bin/bash 

cd bams 
cp ../*.pl .

KK=`ls *bam | grep -v ess | grep -v rmdup | grep -v filt | grep -v deseq | sed "s/.bam//g"` ## list of samples 

for i in $KK 
do
  echo "Samtools: indexing file $i.."
  samtools index $i.bam &
done 
wait 
echo "Done indexing BAM files!" 
echo

for i in $KK
do
  while [ $(jobs | wc -l) -ge 4 ] ; do sleep 5; done
  echo "Picard: marking duplicates in paired-end aligment $i.."
  ## no parallel processing here since Picard is extra annoying on our cluster
  picard MarkDuplicates I=$i.bam O=$i.rmdup.bam M=$i.metrics.out REMOVE_DUPLICATES=true &> $i.rmdup.log & 
done 
wait
echo "Done marking duplicates in BAM files!"
echo

for i in $KK
do
  while [ $(jobs | wc -l) -ge 8 ] ; do sleep 5; done
  echo "Making filtered BAM: processing file $i.."
  ## - samtools only keeps aligned R1 reads and discards unaligned ; 
  ## - cigar_filter.pl removes all reads that don't have salmonella seq right after the tag; 
  samtools view -@4 -h -f 64 -F 4 $i.rmdup.bam | ./cigar_filter.pl | samtools view -@4 -b - > $i.filt.bam & 
done 
wait
echo "Done making filtered BAM files!"  
echo

for i in $KK 
do
  while [ $(jobs | wc -l) -ge 8 ] ; do sleep 5; done
  echo "Making DESeq BAM: processing file $i.."
  ## - take filtered BAM, convert each read to 1nt   
  samtools view -@4 -h $i.filt.bam | ./make_1nt.pl | samtools view -@4 -b - > $i.deseq.bam & 
done 
wait 
echo "Done making DESeq BAM files!"  
echo

## this section is done assuming you have multiple replicates of each condition 
## (here names are *_1.bam, *_2.bam, etc) 
## for essentiality analysis, we merge replicates to gain power

PP=`echo $KK | tr ' ' '\n' | perl -ne 's/_\d+$//g; print' | sort | uniq`

for j in $PP
do
  while [ $(jobs | wc -l) -ge 8 ] ; do sleep 5; done
  echo "Merging BAM files for essentiality analysis: making merged file for $j.."
  samtools merge -@4 $j.merged.bam $j*.filt.bam & 
done 
wait 
echo "Done making merged BAM files!"  
echo

for j in $PP 
do 
  while [ $(jobs | wc -l) -ge 8 ] ; do sleep 5; done
  echo "Making essentiality BAM files for sample $j.."
  samtools view -@4 -h $j.merged.bam | ./make_1nt_uniq.pl | samtools view -@4 -b - > $j.ess.bam & 
done
wait 
echo "Done making merged essentiality BAM files!"
echo

mv *out *log ../logs/
rm *pl
 
echo "ALL BAM FILES WRITTEN SUCCESSFULLY!" 
