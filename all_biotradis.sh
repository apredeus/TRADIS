#!/bin/bash

## this is to get essentiality scores

EMBL=$1
FA=$2
STR=${EMBL%%_*}
EMBL=`readlink -f $EMBL`
FA=`readlink -f $FA`

source activate tradis 

cd fastqs 
KK=`ls *.R1.fastq.gz | grep $STR`

for i in $KK
do
  echo "Adding fake tag TAAGAGACAG; processing sample $i"
  zcat $i | awk '{if (NR%4==2) {print "TAAGAGACAG"$0} else if (NR%4==0) {print "AAAAAAAAAA"$0} else {print}}' | gzip - > ${i%%.R1.fastq.gz}.tag.fastq.gz &
done
wait 

for i in $STR*.tag.fastq.gz 
do
  echo $i > $i.txt
done

for i in $STR*.tag.fastq.gz
do
  echo "Running bacteria_tradis pipeline; processing sample $i" 
  bacteria_tradis -n 16 -v --smalt_r 0 -m 0 -f $i.txt -t TAAGAGACAG -r $FA &> $i.log &
done 
wait 

mv $STR*.tag.* ../biotradis
cd ../biotradis

PP=`ls *insert_site_plot.gz | grep $STR`
for i in $PP
do
  echo "Processing insertion file $i"
  tradis_gene_insert_sites -trim3 0.1 $EMBL $i
done

echo "ALL DONE!" 
