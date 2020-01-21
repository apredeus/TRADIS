# TraDIS processing pipeline
In-house pipeline for processing of TIS/TRADIS data. Version 2 has Bio::Tradis dependency removed, and does not require annotation in EMBL format. 

These scripts allow to perform essentiality analysis, relative fitness analysis, and visualization in JBrowse. Based on the custom pipeline developed by Roy Chaudhuri (r.chaudhuri at sheffield.ac.uk).

## Author
[Alexander Predeus](https://www.researchgate.net/profile/Alexander_Predeus), [Jay Hinton Laboratory](http://www.hintonlab.com/), [University of Liverpool](https://www.liverpool.ac.uk/)

(c) 2019-2020, GPL v3 license

## Dependency installation and getting started 

To get started, you need to install: bwa, featureCounts/Subread, cutadapt, samtools, bedtools, Picard tools, and deepTools. You will also need Perl and R with DESeq2, ggplot2, MASS, and crayon libraries installed. It is convenient to install most of the terminal dependencies using bioconda (R and appropriate libraries might need to be installed separately).

```bash 
conda create -n tradis
source activate tradis
conda install bwa
conda install cutadapt
conda install subread  
conda install samtools 
conda install bedtools
conda install picard 
conda install deeptools 
```

Choose a working directory with enough space (WDIR) and run 
```bash 
cd $WDIR
mkdir raw fastqs bams counts stats browser logs 
```

Put your HiSeq data (e.g. the whole lane of HiSeq sequencing) into WDIR/raw. Clone this repo and copy all the shell/Perl scripts into WDIR. Get your strain (or strains) reference genome and annotation in GFF3 format and place it in your WDIR. Get your barcodes and make a fasta file with them, name it **barcodes.fa**. Your barcodes should contain both Illumina index and transposon sequence, as illustrated below: 

<img src="https://github.com/apredeus/TRADIS/blob/master/img/barcodes.png">

Hence, barcodes file should look like this (each sequence named after the sample it is found in): 

~~~barcodes.fa
>P125109_input_LB_1
TACAAGGCTTCAGGGTTGAGATGTGTATAAGAGACAG
>P125109_input_LB_2
GATCTGGCTTCAGGGTTGAGATGTGTATAAGAGACAG
>P125109_output_LB_10h_1
CGTGATGCTTCAGGGTTGAGATGTGTATAAGAGACAG
>P125109_output_LB_10h_2
CACTGTGCTTCAGGGTTGAGATGTGTATAAGAGACAG
>P125109_output_MAC_1
AAGCTAGCTTCAGGGTTGAGATGTGTATAAGAGACAG
>P125109_output_MAC_2
TGGTCAGCTTCAGGGTTGAGATGTGTATAAGAGACAG
~~~

## Reference and annotation

Experiments for several strains can be processed at once. Each strain needs to have genome assembly (STR.fa) and annotation in GFF3 format (STR.gff). Each meaningful feature should have a locus tag. [Prokka GFF](https://github.com/tseemann/prokka/) would work perfectly well. 

Annotation in GFF3 format needs to be processed in such way that all the features of interest are listed as a **gene**. If your original GFF file was called STR.gff, the new file will be named STR.tradis.gff. Additionally, we will make another GFF file for essentiality analysis - same as STR.tradis.gff, but with last 10% of the gene removed (equivalent of -trim3 0.1 options in Bio-Tradis). This file would be named STR.ess_90.gff. This (and bwa index building) is done by running the following script: 

`./make_refs.sh`

## Demultiplexing

Demultiplexing is unfortunately not parallelized, so it's quite slow. Run it like this, no arguments: 

`./run_cutadapt.sh` 

Cutadapt options used are as follows. Adapter is expected on the 5' (-g), and at least 34 nt have to be matched (otherwise you get garbage matches, e.g. 3 nt long). Default error rate is 0.1, so 3 mismatches over 37 nt sequence are allowed. 

```bash
cutadapt -O 34 -g file:$BC $R1 $R2 -o fastqs/$TAG.R1.fastq.gz -p fastqs/$TAG.R2.fastq.gz --discard-untrimmed
```

This should generate the demultiplexed fastq files in WDIR/fastqs.

## Alignment of demultiplexed reads

Running make_refs.sh should have already built bwa index of all reference genomes (STR.fa). If not, create them by running

`bwa index STR.fa`

After this, alignment is done using **bwa**. Run it with 

`./run_bwa.sh STR.fa`

Actual commands being run are as follows: 

```bash 
bwa mem -t 16 $REF $R1 $R2 | samtools sort -O BAM - > $TAG.bam
```
This generates one sorted paired-end BAM file per sample. More processing will follow.

## Additional BAM processing 

As a part of TRADIS processing, following steps are taken after the alignment:

* Picard tools are used to remove PCR/optical duplicates, generating *Sample.rmdup.bam*;
* Samtools used to keep only aligned R1 reads, and *cigar_filter.pl* used to remove reads that don't have Salmonella sequence right after the TRADIS tag (ie if first bp is not aligned). This generates files named *Sample.filt.bam*; 
* Script named *make_1nt.pl* is used to make 1nt reads located **exactly at the transposon insertion site**. This should minimize quantification ambiguity, and allow to make individual insertion track for visualization. These files are called *Sample.deseq.bam*, and are used for subsequent DESeq2 analysis;
* Script named *make_1nt_uniq.pl* is used to select 1nt reads and keep 1 such read per each insertion site. These files are called *Sample.ess.bam*, and are used for essentiality analysis. 

All these steps should be performed when the following script is ran: 
`./make_bam_files.sh`

## Read counting

One-nucleotide single-end read BAM file (*Sample.deseq.bam*) are used to get gene counts used in DESeq2 (fitness) analysis, while *Sample.ess.bam* is used for essentiality analysis. Two reference GFF files (STR.tradis.gff, STR.ess_90.gff) prepared as shown above, are used for fitness and essentiality analysis, respectively. Quantification is not strand-specific, is calculated per gene feature, and summarized per locus tag (ID). Make sure the logs say that the number of features and meta-features is the same in your GFF, and equals to the number of genes you expect to get. Multimapping reads and reads mapping to more than 1 feature are counted either as-is (for essentiality analysis), or as a fraction (for fitness analysis). 

Actual commands being run: 

```bash
featureCounts -M -O --fraction -t gene -g ID -s 0 -a STR.tradis.gff -o Sample.deseq.fc.tsv Sample.deseq.bam
featureCounts -M -O            -t gene -g ID -s 0 -a STR.ess_90.gff -o Sample.ess.fc.tsv   Sample.ess.bam
```
To get all gene counts, simply run (for each strain): 

`./run_fcount.sh STR`

## Essentiality and fitness analysis 

This step of analysis is done in R. Make sure you have all the needed libraries installed. Copy all results. Make tab-separated table annotation file for each strain, and conditions file (see Conditions.txt in /example). 

*example_analysis.R* gives an example of DESeq2 fitness analysis and essentiality analysis. 

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

## Coverage tracks for visualization

In order to get visualization, you need coverage tracks in bigWig format. To obtain these, run 

`./make_bigwigs.sh`

Actual commands being run: 

```bash
bamCoverage -b Sample.filt.bam -bs 1 -o Sample.bw 
```
## Unique insertion sites and essentiality GFFs

Visualization of coverage, essentiality, and unique insertion sites is conveniently done in JBrowse. Two latter tracks are best presented as GFF files. 

In order to generate these, run the following commands: 

```bash 
./make_ins_gff.sh 
./make_ess_gff.sh STR 
```
Resulting files should generate the following picture once properly set up in JBrowse:

<img src="https://github.com/apredeus/TRADIS/blob/master/img/jbrowse.png">
