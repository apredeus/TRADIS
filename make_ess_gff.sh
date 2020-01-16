#!/bin/bash 

## this one makes GFF files with genes classified by essentiality 

STR=$1        ## eg D7795
GFF=$2        ## eg D7795_v3_2018.tradis.gff 

cd bams
KK=`ls *.ess.bam | grep $STR | sed "s/.ess.bam//"`
cd ..

for i in $KK
do
  ./essential_gff.pl $i $STR.ann_ess.tsv $GFF > $i.ess.gff
done 

mv *.ess.gff browser

echo "MAKING ESSENTIALITY GFF FILES: ALL DONE!" 
