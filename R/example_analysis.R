stringsAsFactors = F
library("MASS")
library("DESeq2")
library("ggplot2")
library("crayon")

## DESeq2 the analysis for D7795 strain 
source("~/TRADIS/R/tradis_functions.R")
setwd("~/TRADIS/example")

exp    <- make_deseq_table("D7795")                    ## read tsv files obtained by featureCounts
cond   <- read.delim("Conditions.txt",row.names=1)     ## read condition file
ann    <- read.delim("D7795_v3_2018.ann",row.names=1)  ## read table-formatted annotation - generated from GFF3 using make_table_ann.pl 

dds    <- DESeqDataSetFromMatrix(countData=round(exp,0), colData=cond, design = ~ Condition)
deseq  <- DESeq(dds)
res1   <- results(deseq, contrast=c("Condition","Output","Input"))
res2   <- results(deseq, contrast=c("Condition","Macro","Input"))
res3   <- results(deseq, contrast=c("Condition","Macro","Output"))
pfc1   <- as.data.frame(res1)[,c(2,6)]
pfc2   <- as.data.frame(res2)[,c(2,6)]
pfc3   <- as.data.frame(res3)[,c(2,6)]
colnames(pfc1)  <- c("input_vs_output_log2FC","input_vs_output_adjP")
colnames(pfc2)  <- c("input_vs_macro_log2FC","input_vs_macro_adjP")
colnames(pfc3)  <- c("output_vs_macro_log2FC","output_vs_macro_adjP")

deseq_table <- merge(ann,exp,by="row.names")
rownames(deseq_table) <- deseq_table$Row.names
deseq_table$Row.names <- NULL
deseq_table <- merge(deseq_table,pfc1,by="row.names")
rownames(deseq_table) <- deseq_table$Row.names
deseq_table$Row.names <- NULL
deseq_table <- merge(deseq_table,pfc2,by="row.names")
rownames(deseq_table) <- deseq_table$Row.names
deseq_table$Row.names <- NULL
deseq_table <- merge(deseq_table,pfc3,by="row.names")

write.table(deseq_table,file="D7795.ann_deseq.tsv",sep="\t",quote=F,row.names=F)

## essenitality analysis - see calculate_essentiality function for details 

ess_table <- make_ess_table("D7795",200)
ess_table <- merge(ann[,c(1,6)],ess_table,by="row.names")
write.table(ess_table,file="D7795.ann_ess.tsv",sep="\t",quote=F,row.names=F)
