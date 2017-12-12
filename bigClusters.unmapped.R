#script to convert results of CD HIT to an R-workspace
#Anna Heintz-Buschart, October 2017

a <- read.delim("validUnmappedSequences.100.100.clustered.bak.clstr",stringsAsFactors=F,header=F)
ar <- a[!grepl("%$",a$V2),]
ar$name <- gsub("... *","",gsub(".+>","",ar$V2),fixed=T)
cc <- read.delim("validUnmappedSequences.100.100.clustered.100.oh.bak.clstr",stringsAsFactors=F,header=F)
cc$name <- gsub("[.]{3} .+","",gsub(".+>","",cc$V2))
arc <- merge(ar[,-2],cc[,-2],by="name",all=T,suffixes=c(".100on100",".bigger"))
a$sample <- gsub(".+>","",gsub("_trim.+","",a$V2))
a$readcount <- as.numeric(gsub("[[:punct:]]{3} [[:punct:]]$","",gsub("... at.+","",gsub(".+-","",a$V2))))
b <- tapply(a$readcount,list(a$V1,a$sample),sum)
b[is.na(b)] <- 0
bd <- merge(b,arc,by.x=0,by.y="V1.100on100",all.x=T)
write.table(bd,"unmappedReadCountsPerCluster100.100PerSample_withBiggerCluster.tsv",sep="\t",quote=F)
saveRDS(bd,"unmappedReadCountsPerCluster100.100PerSample_withBiggerCluster.RDS")
