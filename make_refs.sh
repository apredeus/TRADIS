#!/bin/bash 

for i in *.gff
do
  TAG=${i%%.gff}
  echo "Reference gff file $i: generating tradis and essentiality GFF files.."
  perl -ne '@t=split/\t/; s/\t$t[2]\t/\tgene\t/; print' $TAG.gff > $TAG.tradis.gff
  perl -ne '@t=split/\t/; $l=$t[4]-$t[3]+1; $off = sprintf "%d",0.1*$l; if ($t[6] eq "+") {$end=$t[4]-$off; s/\t$t[4]\t/\t$end\t/; print} else {$beg=$t[3]+$off; s/\t$t[3]\t/\t$beg\t/; print}' $TAG.tradis.gff > $TAG.ess_90.gff
done

for i in `ls *fa | grep -iv barcodes`
do
  echo "Reference fasta file $i: generating bwa indexes.."
  bwa index $i
done 

echo "MAKING REFERENCES: ALL DONE!" 
