perl -ne 's/\tCDS\t/\tgene\t/g; s/\tncRNA\t/\tgene\t/g; s/\trRNA\t/\tgene\t/g; s/\ttRNA\t/\tgene\t/g; s/\ttmRNA\t/\tgene\t/g; print' D7795_v1_2018.gff > D7795_v1_2018.tradis.gff 
perl -ne 's/\tCDS\t/\tgene\t/g; s/\tncRNA\t/\tgene\t/g; s/\trRNA\t/\tgene\t/g; s/\ttRNA\t/\tgene\t/g; s/\ttmRNA\t/\tgene\t/g; print' P125109_v1_2018.gff > P125109_v1_2018.tradis.gff

cat D7795_v1_2018.tradis.gff | sed "s/ID=//g" | sed "s/;Name=/\t/g" | sed "s/;Dbxref=/\t/g" | sed "s/;note=/\t/g" | sed "s/;gene_biotype=/\t/g" | awk -F "\t" '{print $9"\t"$10"\t"$1"\t"$4"\t"$5"\t"$7"\t"$12}' > D7795_v1_2018.ann.tsv
cat P125109_v1_2018.tradis.gff | sed "s/ID=//g" | sed "s/;Name=/\t/g" | sed "s/;Dbxref=/\t/g" | sed "s/;note=/\t/g" | sed "s/;gene_biotype=/\t/g" | awk -F "\t" '{print $9"\t"$10"\t"$1"\t"$4"\t"$5"\t"$7"\t"$12}' > P125109_v1_2018.ann.tsv
 
 
