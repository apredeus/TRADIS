make_tradis_table = function (dir) {
  tsvs       <- list.files(path=dir, pattern="*.fc.tsv")
  cat(sprintf("Adding expression table 1, file name %s\n",tsvs[1]))
  counts     <- read.table(tsvs[1],sep="\t",header=T)
  row.names(counts) <- counts$Geneid
  counts  <- counts[,7,drop=F]
  for (i in 2:length(tsvs)) {
    cat(sprintf("Adding expression table %s, file name %s\n",i,tsvs[i]))
    tmp   <- read.table(tsvs[i],sep="\t",header=T)
    name  <- colnames(tmp)[7]
    counts[[name]] <- tmp[,7]
  }
  colnames(counts) <- gsub(".1nt_rmdup.bam","",colnames(counts))
  colSums(counts)
  return(counts)
}
