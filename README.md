# TRADIS
In-house pipeline for processing of TIS/TRADIS data

This pipeline is based on Bio::Tradis and DESeq2 and allows to perform essentiality analysis, relative fitness analysis, and visualization in JBrowse. Based on the work of Roy Chaudhuri (r.chaudhuri at sheffield.ac.uk). 

## Author
[Alexander Predeus](https://www.researchgate.net/profile/Alexander_Predeus), [Jay Hinton Laboratory](http://www.hintonlab.com/), [University of Liverpool](https://www.liverpool.ac.uk/)

(c) 2019, GPL v3 license

## Getting started 

To get started, you need to install Bio::Tradis, bwa, cutadapt, samtools, and bedtools. You will also need R/Rscript with DESeq2 installed. 

```bash 
conda create -n tradis python=3.6
source activate tradis
conda install cutadapt
conda install bwa 
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

Get your strain's reference genome and GTF and place it in your WDIR. 

Get your barcodes and make a fasta file with them, name it barcodes.fa. Your barcodes should contain both Illumina index and transposon sequence, as illustrated below: 

<IMG> 

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

Demultiplexing is unfortunately not parallelized, so it's quite slow. Run it with 
`./run_cutadapt.sh` 
Cutadapt options used are as follows. Adapter is expected on the 5' (-g), and at least 34 nt have to be matched (otherwise you get garbage matches of e.g. 3 nt). Default error rate is 0.1, so 3 mismatches over 37 nt sequence are allowed. 

```bash
cutadapt -O 34 -g file:$BC $R1 $R2 -o fastqs/{name}.R1.fastq.gz -p fastqs/{name}.R2.fastq.gz --discard-untrimmed
```

This should generate the demultiplexed fastq files in WDIR/fastqs. 
