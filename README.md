# TRADIS
In-house pipeline for processing of TIS/TRADIS data

This pipeline is based on Bio::Tradis and DESeq2 and allows to perform essentiality analysis, relative fitness analysis, and visualization in JBrowse. Based on the custom pipeline developed by Roy Chaudhuri (r.chaudhuri at sheffield.ac.uk). 

## Author
[Alexander Predeus](https://www.researchgate.net/profile/Alexander_Predeus), [Jay Hinton Laboratory](http://www.hintonlab.com/), [University of Liverpool](https://www.liverpool.ac.uk/)

(c) 2019, GPL v3 license

## Getting started 

To get started, you need to install: 

* Bio::Tradis
* bwa
* featureCounts (subread package)
* cutadapt
* samtools
* bedtools
* Picard tools
* GFF

You will also need Perl and R/Rscript with DESeq2 installed.

It is convenient to install most of the dependencies using bioconda. 

```bash 
conda create -n tradis python=3.6
source activate tradis
conda install cutadapt
conda install bwa
conda install subread  
conda install samtools 
conda install bedtools
conda install picard 
```

Choose a working directory with enough space (WDIR) and run 
```bash 
cd $WDIR
mkdir raw fastqs bams counts biotradis browser logs 
```

Put your HiSeq data (e.g. the whole lane of HiSeq sequencing) into WDIR/raw. Clone this repo and copy all the scripts into WDIR. 

Replace tradis_essentiality.R (find it with `which tradis_essentiality.R`) with one from my repo. 

Get your strain's reference genome and GFF and place it in your WDIR. 

Get your barcodes and make a fasta file with them, name it **barcodes.fa**. Your barcodes should contain both Illumina index and transposon sequence, as illustrated below: 

<img src="https://github.com/apredeus/TRADIS/blob/master/img/barcodes.png">

~~~barcodes.fa
>P125109_input_LB_1
CGTGATGCTTCAGGGTTGAGATGTGTATAAGAGACAG
>P125109_input_LB_2
GCCTAAGCTTCAGGGTTGAGATGTGTATAAGAGACAG
>P125109_output_LB_10h_1
ACATCGGCTTCAGGGTTGAGATGTGTATAAGAGACAG
>P125109_output_LB_10h_2
TGGTCAGCTTCAGGGTTGAGATGTGTATAAGAGACAG
~~~

## Annotation formats 

There are many flavours of GFF annotation files, and they all stink. There is an impressive lack of any standards about what has to be included into GFF file. In order to for this pipeline to run smoothly, you would need two things:

* tradis-friendly version of the GFF file
* same annotation in EMBL format

In order to get the first one, you need a GFF in which 1) all the features you care about will be listed as a "gene"; 2) each gene is identified by its locus_tag, e.g. ID=STMMW_00001, or locus_tag=STMMW_00001. 

The exact conversion command would be dependent on the starting GFF. Example of commands I've used for Enteritidis strains: 

```bash 
perl -ne 's/\tCDS\t/\tgene\t/g; s/\tncRNA\t/\tgene\t/g; s/\trRNA\t/\tgene\t/g; s/\ttRNA\t/\tgene\t/g; s/\ttmRNA\t/\tgene\t/g; print' D7795_v1_2018.gff > D7795_v1_2018.tradis.gff 
perl -ne 's/\tCDS\t/\tgene\t/g; s/\tncRNA\t/\tgene\t/g; s/\trRNA\t/\tgene\t/g; s/\ttRNA\t/\tgene\t/g; s/\ttmRNA\t/\tgene\t/g; print' P125109_v1_2018.gff > P125109_v1_2018.tradis.gff
```

For EMBL-based annotation, you would need a tool named [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3), or any other that works for you. Installation should be as simple as `pip install --user git+https://github.com/NBISweden/EMBLmyGFF3.git`, but don't do it inside *tradis* virtual environment - pip often 'sees' other python dependencies outside of the env, which breaks everything. 

Run the conversion using the following command: 

```bash 
EMBLmyGFF3 STR.tradis.gff STR.fa > STR.embl
```
You will be asked a series of questions. Most of these don't make any difference. Use empty locus_tag (program would add XXX\_ to already existing ones!), and '11' for translation table. After this, replace XXX\_ with empty string using your favourite text editor.

If you want gene common names included in the CSV produced by Bio::Tradis, make sure you have 'gene=thrA' entries in the original GFF. 

## Demultiplexing

Demultiplexing is unfortunately not parallelized, so it's quite slow. Run it like this, no arguments: 

`./run_cutadapt.sh` 

Cutadapt options used are as follows. Adapter is expected on the 5' (-g), and at least 34 nt have to be matched (otherwise you get garbage matches of e.g. 3 nt). Default error rate is 0.1, so 3 mismatches over 37 nt sequence are allowed. 

```bash
cutadapt -O 34 -g file:$BC $R1 $R2 -o fastqs/{name}.R1.fastq.gz -p fastqs/{name}.R2.fastq.gz --discard-untrimmed
```

This should generate the demultiplexed fastq files in WDIR/fastqs.

## Alignment and duplicate removal 

Make bwa index of the reference genome (STR.fa): 

`bwa index STR.fa`

After this, alignment is done using **bwa**, and duplicates are removed using **picard MarkDuplicates**. Run it with 

`./all_bwa.sh STR.fa`

Actual commands being run are as follows: 

```bash 
bwa mem -t 16 $REF $R1 $R2 | samtools sort -O BAM - > $TAG.bam
samtools index $TAG.bam
picard MarkDuplicates I=$TAG.bam O=$TAG.rmdup.bam M=$TAG.metrics.out REMOVE_DUPLICATES=true &> $TAG.rmdup.log
```
This generates two BAM files per each sample - one simple alignment, and one with duplicates removed. More processing will follow.

## Additional BAM processing 

As a part of TRADIS processing, following things are done: 

* pairchifilter.pl removes all reads that don't start right after the TRADIS tag (ie if first bp is not aligned)
* samtools view -f 64 -q 10 removes multimappers and R2 reads
* long perl one-liner converts reads into 1nt sequences so that quantification would be less ambiguous (ie if there's an operon or gene overlap)

Actual commands being run: 

```bash
samtools view -h $i | ./pairchifilter.pl | samtools view -h -f 64 -q 10 - | perl -ne 'if (m/^@/) {print} else {@F = split/\t/; $F[5]="1M"; $r=1 if $F[1]&0x10; $F[3]+=length($F[9])-1 if $r; $F[9]=substr $F[9],$r?0:-1,1; $F[10]=substr $F[10],$r?0:-1,1; print join "\t", @F}' | samtools view -b - > ${i%%.rmdup.bam}.1nt_rmdup.bam
samtools view -h $i | ./pairchifilter.pl | samtools view -h -f 64 -q 10 - | samtools view -b - > ${i%%.rmdup.bam}.filt.bam 
```

## Quantification 

To get quantification, you need to convert your GFF into tradis-friendly format. Basically, make sure all the features you are interested in are called **gene**, and that each of those has an ID=locus_tag record. The re-formatted file should be called **STR.tradis.gff**. 

To get all gene counts, simply run

`./all_fc.sh STR.tradis.gff`

Actual commands being run: 

```bash 
featureCounts -t gene -g ID -s 0 -a $GFF -o ${i%%.1nt_rmdup.bam}.fc.tsv $i
```
One-nucleotide single-end read BAM file (\*.1nt_rmdup.bam) are used to get gene counts. Quantification is not strand-specific, is calculated per gene feature, and summarized per locus tag (ID). Make sure the logs say that the number of features and meta-features is the same in your GFF, and equals to the number of genes you expect to get. 

## Coverage tracks 

In order to get visualization, you need coverage tracks in bigWig format. To obtain these, we run 

`./all_bigwig.sh`

Actual commands being run: 

```bash 
bedtools genomecov -ibam $i -bg > ${i%%.filt.bam}.bedGraph
bedGraphToBigWig $i $CHROM ${i%%.bedGraph}.bw
```

CHROM is the STR.chrom.sizes file (tab separated chrom name - chrom size), generated from your reference fasta, STR.fa. 

## Biotradis

Now to the fun part. In order to get biotradis processing to work, and to be consistent with the rest of the pipeline, the following tricks were used: 

* artificial tag TAAGAGACAG is added to every cutadapt-demultiplexed fastq file 
* bacteria_tradis then ran on the resulting (single-end) files to match the introduced tag and to align the reads
* tradis_gene_insert_sites ran to identify unique insertion sites 
* modified tradis_essentiality.R ran on the resulting CSV files to get the essentiality for each gene

All this should be accomplished by running 

`./all_biotradis.sh STR.embl STR.fa`

Actual commands being run:

```bash 
zcat $i | awk '{if (NR%4==2) {print "TAAGAGACAG"$0} else if (NR%4==0) {print "AAAAAAAAAA"$0} else {print}}' | gzip - > ${i%%.R1.fastq.gz}.tag.fastq.gz
echo ${i%%.R1.fastq.gz}.tag.fastq.gz > ${i%%.R1.fastq.gz}.tag.fastq.gz.txt
bacteria_tradis -n 16 -v --smalt_r 0 -m 0 -f $i.txt -t TAAGAGACAG -r $FA
tradis_gene_insert_sites -trim3 0.1 $EMBL $i
tradis_essentiality.R $i
```
This should generate pdf plots of the essentiality analysis distributions and QC, and whole lot of other files. You need \*.tradis_gene_insert_sites.csv as well as \*.param.out file.

## Making GFF files of unique insertion and essentiality 

Visualization is done using JBrowse, so the easiest way to get good feature annotation is to convert it to GFF file (long BED would do too, but we chose GFF). 

In order to get these, run the following commands: 

```bash 
./make_ins_gff.sh STR.fa 
./make_ess_gff.sh STR.tradis.gff 
```

Actual commands being run: 

```bash 
bedtools bamtobed -i $i | awk '{if ($6=="+") {print $1"\t"$2"\t"$2"\t+"} else {print $1"\t"$3"\t"$3"\t-"}}' | sort -k1,1 -k2,2n | uniq | awk '{print $1"\t"$2"\t"$3"\ttag_"NR"\t.\t"$4}' | perl -ne 'chomp; @t=split/\t/; if ($t[5] eq "+") {printf "%s\tTRADIS\tCDS\t%d\t%d\t.\t+\t.\tID=%s\n",$t[0],$t[1]+1,$t[1]+1,$t[3]} else {printf "%s\tTRADIS\tCDS\t%d\t%d\t.\t-\t.\tID=%s\n",$t[0],$t[1],$t[1],$t[3]}' > ${i%%.filt.bam}.uniq_ins.gff
./essential_gff.pl ${i%%.tag.out.gz*} $GFF > $TAG.essential.gff
```

Resulting files should generate the following picture once properly set up in JBrowse:

<img src="https://github.com/apredeus/TRADIS/blob/master/img/jbrowse.png">

## DESeq2 analysis

These steps are done in R. Place all obtained outputs of featureCounts (located in WDIR/counts/\*fc.tsv) into a data directory. 

In the same directory, place the file **Conditions.txt** that would be used to set up comparisons in DESeq2. File should look like this: 

~~~Conditions.txt
Sample	Condition
D23_Input_LB_1	Input
D23_Input_LB_2	Input
D23_RAWpass1_1	RAW_pass1
D23_RAWpass1_2	RAW_pass1
D23_RAWpass3_1	RAW_pass3
D23_RAWpass3_2	RAW_pass3
~~~

After this is done, use the following code (documented in this repo, R subdir): 

```R
stringsAsFactors = F
setwd("~/path/to/my/data")
source("~/path/to/make_tradis_table.R") 

library("MASS")
library("DESeq2")

exp    <- make_tradis_table(getwd())
cond   <- read.table("Conditions.txt",header=T,row.names=1)

dds    <- DESeqDataSetFromMatrix(countData=round(exp,0),colData=cond,design = ~ Condition)
deseq  <- DESeq(dds)
res    <- results(deseq, contrast=c("Condition","RAW_pass1","Input"))
```

This should generate the table with DESeq2 data. 
