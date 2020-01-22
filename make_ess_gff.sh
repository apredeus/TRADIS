#!/bin/bash 

## this one makes GFF files with genes classified by essentiality 
## you need to run R section of the analysis and obtain $PREFIX.ann_ess.tsv table for each strain

STR=$1        ## eg D7795_v3_2018
PREFIX=${STR%%_*}

cd bams
KK=`ls *.ess.bam | grep $PREFIX | sed "s/.ess.bam//"`
cd ..

for i in $KK
do
  ./essential_gff.pl $i $PREFIX.ann_ess.tsv $STR.tradis.gff > $i.ess.gff
done 

mv *.ess.gff browser

echo "MAKING ESSENTIALITY GFF FILES: ALL DONE!" 
