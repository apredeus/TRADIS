# TRADIS
In-house pipeline for processing of TIS/TRADIS data

This pipeline is based on Bio::Tradis and DESeq2 and allows to perform essentiality analysis, relative fitness analysis, and visualization in JBrowse. Based on the work of Roy Chaudhuri (r.chaudhuri at sheffield.ac.uk). 

## Author
[Alexander Predeus](https://www.researchgate.net/profile/Alexander_Predeus), [Jay Hinton Laboratory](http://www.hintonlab.com/), [University of Liverpool](https://www.liverpool.ac.uk/)

(c) 2019, GPL v3 license

## Getting started 

To get started, you need to install Bio::Tradis, bwa, featureCounts (from subread package), cutadapt, samtools, and bedtools. You will also need R/Rscript with DESeq2 installed. 

```bash 
conda create -n tradis python=3.6
source activate tradis
conda install cutadapt
conda install bwa
conda install subread  
conda install samtools 
conda install bedtools
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

Now to the fun part. In order to get biotradis processing to work, and to be consistent with the rest of the pipeline, the following tricks wre used. 

