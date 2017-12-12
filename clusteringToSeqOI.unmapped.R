#R script to convert clustering results of CD-HIT to an R workspace (also prints to .tsv)

a <- read.delim("validSequences.unmapped.100.clustered.14.seqIO.bak.clstr",stringsAsFactors=F,header=F)
aS <- a[grep("*",a$V2,invert=T,fixed=T),]
aS$sample <- gsub(".+>","",gsub("_trim.+","",aS$V2))
aS$readcount <- as.numeric(gsub("... at.+","",gsub(".+-","",aS$V2)))
b <- tapply(aS$readcount,list(aS$V1,aS$sample),sum)
b[is.na(b)] <- 0
aOI <- a[grep("*",a$V2,fixed=T),] 
aOI$seq <- gsub(".+>","",gsub("... *","",aOI$V2,fixed=T))
aOI <- aOI[order(aOI$V1),]
rownames(b) <- aOI$seq
write.table(b,"readCountsPerSeqOIPerSample.unmapped.14.tsv",sep="\t",quote=F)
saveRDS(b,"readCountsPerSeqOIPerSample.unmapped.14.RDS")
