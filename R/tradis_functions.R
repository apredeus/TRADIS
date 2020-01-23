make_deseq_table = function (strain) {
  tsvs       <- list.files(pattern="*.deseq.fc.tsv")
  str_tsvs   <- Filter(function(x) grepl(strain, x), tsvs)
  cat(sprintf("Adding expression table 1, file name %s\n",str_tsvs[1]))
  deseq_matrix     <- read.table(tsvs[1],sep="\t",header=T)
  row.names(deseq_matrix) <- deseq_matrix$Geneid
  deseq_matrix  <- deseq_matrix[,7,drop=F]
  for (i in 2:length(tsvs)) {
    cat(sprintf("Adding expression table %s, file name %s\n",i,str_tsvs[i]))
    tmp   <- read.table(tsvs[i],sep="\t",header=T)
    name  <- colnames(tmp)[7]
    deseq_matrix[[name]] <- tmp[,7]
  }
  colnames(deseq_matrix) <- gsub(".deseq.bam","",colnames(deseq_matrix))
  colSums(deseq_matrix)
  return(deseq_matrix)
}

make_ess_table = function (strain,len_cutoff=200) {
  ## read all appropriate TSV files made for essentiality analysis
  ## and calculate insertion index for each gene 
  ## (by default, only first 90% of the gene are considered - ~ to -trim3 0.1 in Bio-Tradis)
  
  tsvs           <- list.files(pattern="*.ess.fc.tsv")
  str_tsvs       <- Filter(function(x) grepl(strain, x), tsvs)
  cat(sprintf("Adding expression table 1, file name %s\n",str_tsvs[1]))
  cat(sprintf("--------------------------------------------------------------------------\n\n"))
  ess_matrix     <- read.table(str_tsvs[1],sep="\t",header=T,row.names = 1)
  tag            <- gsub(".ess.bam","",colnames(ess_matrix)[6])
  tag_count      <- paste(tag,"ins_count",sep=".")
  colnames(ess_matrix)[6] <- tag_count
  tag_ii         <- paste(tag,"ins_index",sep=".")
  tag_ess        <- paste(tag,"ess",sep=".")
  ## calculate insertion index 
  ess_matrix[[tag_ii]]  <- ess_matrix[[tag_count]]/ess_matrix$Length
  ii_len                <- ess_matrix[,c(tag_ii,"Length")]
  ## use insertion index to calculate essentiality - see function below
  ess_matrix[[tag_ess]] <- calculate_essentiality(ii_len,len_cutoff)  
  cat(sprintf("==========================================================================\n\n"))

  
  for (i in 2:length(str_tsvs)) {
    cat(sprintf("Adding expression table %s, file name %s\n",i,str_tsvs[i]))
    cat(sprintf("--------------------------------------------------------------------------\n\n"))
    tmp   <- read.table(str_tsvs[i],sep="\t",header=T,row.names=1)
    tag            <- gsub(".ess.bam","",colnames(tmp)[6])
    tag_count      <- paste(tag,"ins_count",sep=".")
    colnames(tmp)[6] <- tag_count
    tag_ii         <- paste(tag,"ins_index",sep=".")
    tag_ess        <- paste(tag,"ess",sep=".")
    ess_matrix[[tag_count]] <- tmp[[tag_count]]
    ess_matrix[[tag_ii]]    <- tmp[[tag_count]]/ess_matrix$Length
    ii_len                  <- ess_matrix[,c(tag_ii,"Length")]

    ## use insertion index to calculate essentiality - see function below
    ess_matrix[[tag_ess]] <- calculate_essentiality(ii_len,len_cutoff)  
    cat(sprintf("==========================================================================\n\n"))
  }
  cat(sprintf("Done making essentiality matrix; unique insertions counts for each condition:\n"))
  colSums(ess_matrix[,grep("ins_count",colnames(ess_matrix))])
  
  return(ess_matrix)
}

calculate_essentiality = function (ii_len,len_cutoff) { 
  library(ggplot2)
  library(crayon)
  library(MASS)
  ## take a vector of insertion indexes, return a vector of essential/ambigous/non-essential
  ## pretty much 100% follows the function from Bio::Tradis
  ## version 2 only fits distributions over genes with Length >= len_cutoff (default 200)
  ## this is important for low-density libraries that don't have enough data for short genes 
  
  #identify second maxima (orig code gave wrong results for our data, so check plots)
  ii_df        <- ii_len
  ii_df$ess    <- "NONE"
  colnames(ii_df)[1] <- "ii"
  which_long   <- which(ii_df$Length >= len_cutoff)
  which_short  <- which(ii_df$Length <  len_cutoff)
  ii_vector    <- ii_df[ii_df$Length >= len_cutoff,]$ii
  cat(sprintf("Using length cutoff of %d; %d genes are selected, %d are too short.\n",len_cutoff,length(which_long),length(which_short)))
  
  ii_df[which_short,]$ess <- "short"  ## all genes shorter than len_cutoff by default deemed non-essential
  
  h        <- hist(ii_vector, breaks = 200, plot=F)
  maxindex <- which.max(h$density[4:length(h$density)])
  maxval   <- h$mids[maxindex+3] 
  
  #find inter-mode minima with loess
  nG       <- length(ii_vector)
  r        <- floor(maxval*2000)
  I        <- ii_vector < r/2000
  h1       <- hist(ii_vector[I], breaks=(0:r/2000), plot=F)
  lo       <- loess(h1$density ~ c(1:length(h1$density))) #loess smothing over density
  m        <- h1$mids[which.min(predict(lo))]
  
  cat(sprintf("Overall number of genes in the vector: %d, zeroes: %d, mean insertion index: %f, stdev: %f\n\n",
              nG,sum(ii_vector==0),mean(ii_vector),sd(ii_vector)))
  cat(sprintf("Estimated second maximum in the insertion index distribution: %f\n",maxval))
  cat(sprintf("Calculated Loess smothing minimum of insertion index: %f\n\n",m))
  
  I1       <- ((ii_vector < m)&(ii_vector >= 0))
  I2       <- ((ii_vector >= m)&(ii_vector < h$mids[max(which(h$counts>5))]))
  
  h        <- hist(ii_vector, breaks="FD", plot=F) 
  f1       <- (sum(I1) + sum(ii_vector == 0))/nG
  f2       <- (sum(I2))/nG
  
  d1       <- fitdistr(ii_vector[I1], "exponential")
  d2       <- fitdistr(ii_vector[I2], "gamma") #fit curves
  cat(sprintf("Fitting exponential function using %d values and gamma function using %d values.\n",sum(I1),sum(I2)))
  #plot histogram and fitted curves 
  
  exp_fit  <- as.data.frame(cbind(0:500/4000,f1*dgamma(0:500/4000, 1, d1$estimate[1])))
  gmm_fit  <- as.data.frame(cbind(0:500/4000,f2*dgamma(0:500/4000, d2$estimate[1], d2$estimate[2])))
  colnames(exp_fit) <- c("val","exp")
  colnames(gmm_fit) <- c("val","gmm")
  
  #calculate log-odds ratios to choose thresholds
  log_ratio <- log((pgamma(1:1000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:1000/10000, 1,d1$e[1], 
                  lower.tail=F)))/(pgamma(1:1000/10000, 1,d1$e[1], lower.tail=F)*(1-pgamma(1:1000/10000, d2$e[1],d2$e[2]))) , base=2)
  lrt_df    <- as.data.frame(cbind(1:1000/10000,log_ratio))
  
  p1 <- ggplot(lrt_df,aes(x=V1,y=log_ratio)) + geom_line(size=2) + xlim(0,maxval) + ylim(-20,20) + 
    geom_hline(yintercept = 2,color="red") + geom_hline(yintercept = -2,color="blue") + 
    labs(title="Log-odds ratio for gamma (non-essential) and exp (essential) distributions")
  print(p1)
  
  lower <- max(which(log_ratio < -2))
  upper <- min(which(log_ratio > 2))
  
  essen <- lower/10000
  ambig <- upper/10000
  
  cat(sprintf("Calculated essentiality cutoff: %f, non-essentiality cutoff: %f\n",essen,ambig))
  
  ## make a nice plot 
  
  p2 <- ggplot(ii_df,aes(x=ii)) + geom_histogram(color="black",fill="light blue",bins=100) + xlim(-1e-3,5*maxval) + 
    geom_line(data=exp_fit,aes(x=val,y=exp),color="red",size=1.5) + geom_line(data=gmm_fit,aes(x=val,y=gmm),color="blue",size=1.5) + 
    geom_vline(xintercept = essen, color="black") + geom_vline(xintercept = ambig, color="grey") + 
    labs(title="Histogram of insertion indexes and fitted essential (red) and non-essential (blue) gene distribution")
  print(p2)
  
  ## make and return a char vector with "essential", "ambiguous", and "non-essential" calls 
  
  calls     <- rep("non-essential",nG)
  ess_index <- which(ii_vector <= essen)
  amb_index <- which(ii_vector > essen & ii_vector <= ambig)
  calls[ess_index] <- "essential"
  calls[amb_index] <- "ambiguous"
  ii_df[which_long,]$ess <- calls
  cat(red("Final statistics: "))
  cat(sprintf("%d essential, %d ambiguous, and %d non-essential genes.\n",
              table(ii_df$ess)[2],table(ii_df$ess)[1],table(ii_df$ess)[3]))
  return(ii_df$ess)
}


