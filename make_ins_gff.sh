#!/bin/bash 

## this script makes GFF files to visualize unique TRADIS insertion sites 

cd bams 

for i in *ess.bam
do
  TAG=${i%%.ess.bam}
  echo "Processing sample $i.."
  bedtools bamtobed -i $i | awk '{if ($6=="+") {print $1"\t"$2"\t"$2"\t+"} else {print $1"\t"$3"\t"$3"\t-"}}' | sort -k1,1 -k2,2n | uniq | awk '{print $1"\t"$2"\t"$3"\ttag_"NR"\t.\t"$4}' | perl -ne 'chomp; @t=split/\t/; if ($t[5] eq "+") {printf "%s\tTRADIS\tCDS\t%d\t%d\t.\t+\t.\tID=%s\n",$t[0],$t[1]+1,$t[1]+1,$t[3]} else {printf "%s\tTRADIS\tCDS\t%d\t%d\t.\t-\t.\tID=%s\n",$t[0],$t[1],$t[1],$t[3]}' > $TAG.uniq_ins.gff & 
done 

wait

mv *gff ../browser 

echo "MAKING INSERTION GFFS: ALL DONE!" 
