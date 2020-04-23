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
  #ess_matrix[[tag_ess]] <- calculate_essentiality(ii_len,len_cutoff)  
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
    #ess_matrix[[tag_ess]] <- calculate_essentiality(ii_len,len_cutoff)  
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
  lo       <- loess(h1$density ~ c(1:length(h1$density)),span = 0.5) #loess smothing over density
  preds    <- as.data.frame(cbind(h1$density,predict(lo)))
  colnames(preds) <- c("hist","loess")
  row.names(preds) <- 1:nrow(preds)
  
  p0 <- ggplot(preds, aes(x = as.numeric(rownames(preds)), y = hist)) + geom_point() + geom_line(aes(y = loess), size = 1)
  print(p0)
  
  m        <- h1$mids[which.min(predict(lo))]
  #if (m < 0.005) m <- 0.005  ## for very saturated libraries
  
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

plot_essentiality_histogram = function (tsv,len_cutoff = 200) { 
  library(ggplot2)
  library(patchwork)
  library(scales)
  ## read all essentiality TSVs and plot them using 3 bin sizes 
  ii_df   <- read.table(tsv,sep="\t",header = T,row.names = 1)
  tag            <- gsub(".ess.bam","",colnames(ii_df)[6])
  colnames(ii_df)[6] <- c("Count")
  ii_df$ii    <- ii_df$Count/ii_df$Length
  which_long     <- which(ii_df$Length >= len_cutoff)
  which_short    <- which(ii_df$Length <  len_cutoff)
  cat(sprintf("Using length cutoff of %d; %d genes are selected, %d are too short.\n",len_cutoff,length(which_long),length(which_short)))
  filt_counts    <- ii_df[which_long,]
  
  p1 <- ggplot(filt_counts,aes(x=ii)) + geom_histogram(fill="steelblue3",bins=70) + xlim(-0.005,0.295) +
          ggtitle(paste(tag,", bins = 70",sep = " ")) 
  p2 <- ggplot(filt_counts,aes(x=ii)) + geom_histogram(fill="steelblue3",bins=200) + xlim(-0.005,0.295) +
          ggtitle(paste(tag,", bins = 200",sep = " "))
  p3 <- ggplot(filt_counts,aes(x=ii)) + geom_histogram(fill="steelblue3",bins=500) + xlim(-0.005,0.295) +
          ggtitle(paste(tag,", bins = 500",sep = " "))
  p4 <- ggplot(filt_counts,aes(x=ii)) + geom_density(bw = 0.001) + xlim(-0.005,0.1)
  suppressWarnings(print((p1 / p2 / p3) | p4))
  }

custom_essentiality = function (tsv,m,len_cutoff = 200,simple = F) {
  ## simplified function for shallow libraries
  ## if simple = T, only genes with 0 ii are deemed essential, and the rest are non-essential (short still applies)
  ## otherwise, take m from the density plot in plot_essentiality_histograms and use it as the divider
  library(ggplot2)
  library(crayon)
  library(MASS)
  
  if (! simple) { 
    ii_df          <- read.table(tsv,sep="\t",header = T,row.names = 1)
    tag            <- gsub(".ess.bam","",colnames(ii_df)[6])
    tag_count      <- paste(tag,"ins_count",sep=".")
    tag_ii         <- paste(tag,"ins_index",sep=".")
    tag_ess        <- paste(tag,"ess",sep=".")
    
    colnames(ii_df)[6] <- tag_count
    ii_df[[tag_ii]]    <- ii_df[[tag_count]]/ii_df$Length
    ii_df[[tag_ess]]   <- "short"

    which_long   <- which(ii_df$Length >= len_cutoff)
    which_short  <- which(ii_df$Length <  len_cutoff)
    ii_filt      <- ii_df[which_long,]
    ii_vector    <- ii_filt[[tag_ii]]
    cat(sprintf("Using length cutoff of %d; %d genes are selected, %d are too short.\n",len_cutoff,length(which_long),length(which_short)))
    
    ## no Lo-ass fitting - just use the estimate from the density. For real, it's way better
    ## feed the estimate as "m" variable
    cat(sprintf("Overall number of genes in the vector: %d, zeroes: %d, mean insertion index: %f, stdev: %f\n\n",
                length(ii_vector),sum(ii_vector == 0),mean(ii_vector),sd(ii_vector)))

    nG       <- length(ii_vector)
    h        <- hist(ii_vector,breaks = 200,plot = F) 
    
    I1       <- (ii_vector < m)
    I2       <- ((ii_vector >= m) & (ii_vector < h$mids[max(which(h$counts>5))])) ## don't ask
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
    log_ratio <- log((pgamma(1:1000/10000,d2$e[1],d2$e[2])*(1 - pgamma(1:1000/10000,1,d1$e[1], 
                    lower.tail = F)))/(pgamma(1:1000/10000,1,d1$e[1], lower.tail = F)*(1 - pgamma(1:1000/10000,d2$e[1],d2$e[2]))),base = 2)
    lrt_df    <- as.data.frame(cbind(1:1000/10000,log_ratio))

    lower <- max(which(log_ratio < -2))
    upper <- min(which(log_ratio > 2))
    
    essen <- lower/10000
    ambig <- upper/10000
    
    cat(sprintf("Calculated essentiality cutoff: %f, non-essentiality cutoff: %f\n",essen,ambig))
    
    ## make a nice plot 
    
    pfit <- ggplot(ii_filt,aes_string(x = tag_ii)) + geom_histogram(color = "black",fill = "light blue",bins = 100) + xlim(-1e-3,0.1) + 
      geom_line(data = exp_fit,aes(x = val,y = exp),color = "red",size = 1.5) + geom_line(data = gmm_fit,aes(x = val,y = gmm),color = "blue",size = 1.5) + 
      geom_vline(xintercept = essen,color = "black",size = 1) + geom_vline(xintercept = ambig,color = "grey",size = 1) + 
      ggtitle(paste("Sample: ",tag,"; histogram of insertion indexes and fitted essential (red) and non-essential (blue) gene distribution",sep = ""))
    suppressWarnings(print(pfit))
    
    ## make and return a char vector with "essential", "ambiguous", and "non-essential" calls 
    
    calls     <- rep("non-essential",nG)
    ess_index <- which(ii_vector <= essen)
    amb_index <- which(ii_vector > essen & ii_vector <= ambig)
    calls[ess_index] <- "essential"
    calls[amb_index] <- "ambiguous"
    ii_df[which_long,][[tag_ess]] <- calls
    cat(red("Final statistics: "))
    cat(sprintf("%d short, %d essential, %d ambiguous, and %d non-essential genes.\n",
                table(ii_df[[tag_ess]])[4],table(ii_df[[tag_ess]])[2],table(ii_df[[tag_ess]])[1],table(ii_df[[tag_ess]])[3]))
  } else { 
    ## "simple" case - ignore m; short genes are short, zeros are essential, rest is non-essential
    ii_df          <- read.table(tsv,sep="\t",header = T,row.names = 1)
    tag            <- gsub(".ess.bam","",colnames(ii_df)[6])
    tag_count      <- paste(tag,"ins_count",sep=".")
    tag_ii         <- paste(tag,"ins_index",sep=".")
    tag_ess        <- paste(tag,"ess",sep=".")
    
    colnames(ii_df)[6] <- tag_count
    ii_df[[tag_ii]]    <- ii_df[[tag_count]]/ii_df$Length
    ii_df[[tag_ess]]   <- "short"
    
    which_long   <- which(ii_df$Length >= len_cutoff)
    which_short  <- which(ii_df$Length <  len_cutoff)
    ii_filt      <- ii_df[which_long,]
    ii_vector    <- ii_filt[[tag_ii]]
    nG           <- length(ii_vector)
    
    cat(sprintf("Using length cutoff of %d; %d genes are selected, %d are too short.\n",len_cutoff,length(which_long),length(which_short)))
    cat(red("WARNING: executing simplified version of TRADIS analysis!\n"))
    cat(sprintf("No fitting will be done, only genes with insertion index of 0 are deemed essential.\n"))
    cat(sprintf("Using length cutoff of %d; %d genes are selected, %d are too short.\n",len_cutoff,length(which_long),length(which_short)))
    
  ## feed the estimate as "m" variable
    cat(sprintf("Overall number of genes in the vector: %d, zeroes: %d, mean insertion index: %f, stdev: %f\n\n",
                length(ii_vector),sum(ii_vector == 0),mean(ii_vector),sd(ii_vector)))
   
    ## make and return a char vector with "essential", "ambiguous", and "non-essential" calls 
    
    calls     <- rep("non-essential",nG)
    ess_index <- which(ii_vector == 0)
    calls[ess_index] <- "essential"
    ii_df[which_long,][[tag_ess]] <- calls
    cat(red("Final statistics: "))
    cat(sprintf("%d short, %d essential, 0 ambiguous, and %d non-essential genes.\n",
                table(ii_df[[tag_ess]])[3],table(ii_df[[tag_ess]])[1],table(ii_df[[tag_ess]])[2]))
  }
  return(ii_df)
}

make_custom_table = function (ann,ess_tables) {
  ## annotation and list
  cat(sprintf("\nAdding essentiality table # 1...\n"))
  out_table <- merge(ann[,c(1,6)],ess_tables[[1]],by = "row.names")
  rownames(out_table) <- out_table$Row.names
  out_table$Row.names <- NULL
  
  for (i in 2:length(ess_tables)) {
    cat(sprintf("Adding essentiality table # %d...\n",i))
    out_table <- merge(out_table,ess_tables[[i]][,c(6,7,8)],by="row.names")
    rownames(out_table) <- out_table$Row.names
    out_table$Row.names <- NULL
  }
  cat(sprintf("\nDone making essentiality matrix; unique insertions counts for each condition:\n\n"))
  print(colSums(out_table[,grep("ins_count",colnames(out_table))]))
  
  return(out_table)
}