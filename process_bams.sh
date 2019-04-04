#!/bin/bash 

cd bams 
cp ../pairchifilter.pl .

for i in *.rmdup.bam 
do
  echo "Processing sample $i (1nt BAM files).."
  ## here's what's happening here: 
  ## - pairchifilter.pl removes all reads that don't start right after the TRADIS tag (ie if first bp is not aligned) 
  ## - samtools view -f 64 -q 10 removes multimappers and R2 reads
  ## - long perl one-liner converts reads into 1nt sequences so that quantification would be less ambiguous (ie if there's an operon or gene overlap) 

  samtools view -h $i | ./pairchifilter.pl | samtools view -h -f 64 -q 10 - | perl -ne 'if (m/^@/) {print} else {@F = split/\t/; $F[5]="1M"; $r=1 if $F[1]&0x10; $F[3]+=length($F[9])-1 if $r; $F[9]=substr $F[9],$r?0:-1,1; $F[10]=substr $F[10],$r?0:-1,1; print join "\t", @F}' | samtools view -b - > ${i%%.rmdup.bam}.1nt_rmdup.bam &  
done 
wait 

for i in *.rmdup.bam
do
  echo "Processing sample $i (filtered BAM files)"
  ## everything is the same, but we keep the whole read to get the visualized coverage
  samtools view -h $i | ./pairchifilter.pl | samtools view -h -f 64 -q 10 - | samtools view -b - > ${i%%.rmdup.bam}.filt.bam &  
done
wait 


echo "BAM PROCESSING: ALL DONE!" 
