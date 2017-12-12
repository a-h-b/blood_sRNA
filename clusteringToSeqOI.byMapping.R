#R script to convert clustering results of CD-HIT, combined with data from a sam file, to an R workspace (also prints to .tsv)

a <- read.delim("validSequences.100.clustered.14.seqIO.bak.clstr",stringsAsFactors=F,header=F)
aS <- a[grep("*",a$V2,invert=T,fixed=T),]
aS$sample <- gsub(".+>","",gsub("_trim.+","",aS$V2))
aS$readcount <- as.numeric(gsub("... at.+","",gsub(".+-","",aS$V2)))
aS$readtype <- ""
aS$read <- gsub("[[:punct:][:alpha:]]+","",gsub("-.+","",gsub(".+_","",aS$V2)))
aOI <- a[grep("*",a$V2,fixed=T),]
aOI$seq <- gsub(".+>","",gsub("... *","",aOI$V2,fixed=T))
aOI <- aOI[order(aOI$V1),]
for(s in unique(aS$sample)){
  if(s !="empty41"){ #not done for the empty library without data
    samf <- list.files(path=paste0("../analysis/",s,"/aln_hg38/"),pattern="sam$")
    sam <- read.delim(paste0("../analysis/",s,"/aln_hg38/",samf),header=F,comment.char="@",stringsAsFactors=F,flush=T)[,c(1,2,15)]
    sam$mm <- gsub(".+:","",sam$V15)
    aS$readtype[aS$sample==s&aS$read %in% unique(gsub("-.+","",sam$V1[sam$mm==""]))] <- "nonmapping"
    aS$readtype[aS$sample==s&aS$read %in% unique(gsub("-.+","",sam$V1[sam$mm=="2"]))] <- "2mm"
    aS$readtype[aS$sample==s&aS$read %in% unique(gsub("-.+","",sam$V1[sam$mm=="1"]))] <- "1mm"
    aS$readtype[aS$sample==s&aS$read %in% unique(gsub("-.+","",sam$V1[sam$mm=="0"]))] <- "perfect"
    if(grepl("-",s)){
      samf <- list.files(path=paste0("../analysis/",s,"/aln_salmonella/"),pattern="sam$")
      sam <- read.delim(paste0("../analysis/",s,"/aln_salmonella/",samf),header=F,comment.char="@",stringsAsFactors=F,flush=T)[,c(1,2,15)]
      sam$mm <- gsub(".+:","",sam$V15)
      aS$readtype[aS$sample==s&aS$readtype==""&aS$read %in% unique(gsub("-.+","",sam$V1[sam$mm=="2"]))] <- "2mmSalmonella"
      aS$readtype[aS$sample==s&aS$readtype==""&aS$read %in% unique(gsub("-.+","",sam$V1[sam$mm=="1"]))] <- "1mmSalmonella"
      aS$readtype[aS$sample==s&aS$readtype==""&aS$read %in% unique(gsub("-.+","",sam$V1[sam$mm=="0"]))] <- "perfectSalmonella"
    }
  }}
for(type in c("nonmapping","2mm","1mm","perfect","2mmSalmonella","1mmSalmonella","perfectSalmonella")){
  b <- tapply(aS$readcount[aS$readtype==type],list(aS$V1[aS$readtype==type],aS$sample[aS$readtype==type]),sum)
  b[is.na(b)] <- 0
  rownames(b) <- aOI$seq[aOI$V1 %in% rownames(b)]
  write.table(b,paste0("readCountsPerSeqOIPerSample.",type,".14.tsv"),sep="\t",quote=F)
  saveRDS(b,paste0("readCountsPerSeqOIPerSample.",type,".14.RDS"))
  rm(b)
}
