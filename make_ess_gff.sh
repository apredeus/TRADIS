#!/bin/bash 

## this one makes GFF files with genes classified by essentiality 

GFF=$1
STR=${GFF%%_*}
GFF=`readlink -f $GFF`

cd biotradis
cp ../essential_gff.pl .
## note that you have to replace tradis_essentiality.R with one in the repository to get .param.out files as well 
## or you can just add two lines to it that print it out 

KK=`ls *.param.out | grep $STR`

for i in $KK
do
  TAG=${i%%.tag.out.gz*}
  ./essential_gff.pl $i $GFF > $TAG.essential.gff
done 

mv *.essential.gff ../browser

echo "ALL DONE!" 
