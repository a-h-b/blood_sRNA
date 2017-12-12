a <- read.delim("validUnmappedSequences.100.100.clustered.bak.clstr",stringsAsFactors=F,header=F)
a$sample <- gsub(".+>","",gsub("_hg38.+","",a$V2))
a$readcount <- as.numeric(gsub("[[:punct:]]{3} [[:punct:]]$","",gsub("... at.+","",gsub(".+-","",a$V2))))
b <- tapply(a$readcount,list(a$V1,a$sample),sum)
b[is.na(b)] <- 0
write.table(b,"unmappedReadCountsPerCluster100.100PerSample.tsv",sep="\t",quote=F)
saveRDS(b,"unmappedReadCountsPerCluster100.100PerSample.RDS")
