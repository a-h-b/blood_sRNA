source("anna3tab2venn.R")

#read counts per clusters and membership to confirmed contaminants
clustB <- readRDS("unmappedReadCountsPerCluster100.100PerSample_withBiggerCluster.RDS")
clusStop <- read.delim("stops.sorted.txt",stringsAsFactors = F, header=F)$V1 # Illumina stop primer
clusAlg <- read.delim("algs.sorted.txt",stringsAsFactors = F, header=F)$V1 # sRNA 1
clus8182 <- read.delim("8182s.sorted.txt",stringsAsFactors = F, header=F)$V1 # sRNA 3
clus3869 <- read.delim("3869s.sorted.txt",stringsAsFactors = F, header=F)$V1 # sRNA 5
clus3378 <- read.delim("3378s.sorted.txt",stringsAsFactors = F, header=F)$V1 # sRNA 6
clus8030 <- read.delim("8030s.sorted.txt",stringsAsFactors = F, header=F)$V1 # sRNA 2
clus6846 <- read.delim("6848s.sorted.txt",stringsAsFactors = F, header=F)$V1 # sRNA 4

#keep only datasets of interest
clustB2 <- clustB[!rownames(clustB) %in% clusStop,]
clustBc <- clustB2[,c(2:3,12:15,24:31,40:47,56:65)]
clustBc <- clustBc[rowSums(clustBc[,1:30])>0,]
clustBc1 <- clustBc[,c(3,4,7,8,11,12,15,16,19,20,23,24,27,28,31:32)]
colnames(clustBc1)[1:14] <- gsub("5.+","",colnames(clustBc1)[1:14])
clustBc1 <- clustBc1[,c(3:4,1:2,5:6,9:10,7:8,11:12,13:16)]

#sum reads by large clusters
clustBc1c <- aggregate(clustBc1[!rownames(clustBc1) %in% c(clusAlg,clus8182,clus3378,clus3869,clus8030,clus6846),1:14],
                       list(clustBc1$V1.bigger[!rownames(clustBc1) %in% c(clusAlg,clus8182,clus3378,clus3869,clus8030,clus6846)]),sum)
rownames(clustBc1c) <- clustBc1c[,1]
clustBc1c <- as.matrix(clustBc1c[,-1])

#tables with potential contaminants that are replicated
consTab <- data.frame("B"=apply(clustBc1c[,1:2],1,function(x) all(x>10))&rowSums(clustBc1c[,13:14])<2,
                      "A"=apply(clustBc1c[,3:4],1,function(x) all(x>10))&rowSums(clustBc1c[,13:14])<2,
                      "C"=apply(clustBc1c[,5:6],1,function(x) all(x>10))&rowSums(clustBc1c[,13:14])<2,
                      "E"=apply(clustBc1c[,7:8],1,function(x) all(x>10))&rowSums(clustBc1c[,13:14])<2,
                      "D"=apply(clustBc1c[,9:10],1,function(x) all(x>10))&rowSums(clustBc1c[,13:14])<2,
                      "F"=apply(clustBc1c[,11:12],1,function(x) all(x>10))&rowSums(clustBc1c[,13:14])<2,
                      stringsAsFactors = F)

cMTs <- data.frame("B"=ifelse(apply(clustBc1c[,1:2],1,function(x) all(x>10))&rowSums(clustBc1c[,13:14])<2,
                              clustBc1c[,1],0),
                   "A"=ifelse(apply(clustBc1c[,3:4],1,function(x) all(x>10))&rowSums(clustBc1c[,13:14])<2,
                              clustBc1c[,3:4],0),
                   "C"=ifelse(apply(clustBc1c[,5:6],1,function(x) all(x>10))&rowSums(clustBc1c[,13:14])<2,
                              clustBc1c[,5:6],0),
                   "E"=ifelse(apply(clustBc1c[,7:8],1,function(x) all(x>10))&rowSums(clustBc1c[,13:14])<2,
                              clustBc1c[,7:8],0),
                   "D"=ifelse(apply(clustBc1c[,9:10],1,function(x) all(x>10))&rowSums(clustBc1c[,13:14])<2,
                              clustBc1c[,9:10],0),
                   "F"=ifelse(apply(clustBc1c[,11:12],1,function(x) all(x>10))&rowSums(clustBc1c[,13:14])<2,
                              clustBc1c[,11:12],0),
                   stringsAsFactors = F)

# normalize to spike-in (from qPCR data)
ex <- (colSums(cMTs)*1e6/c(4017825,4239897,3339719,2372020,4103878,3569672))*40000/c(123.2009856,1032.572254,488.963293,266.4395747,222.9598443,527.2193075)

#plot
pdf("Figure4D.pdf",width=7/2.54,height=2.3/2.54,pointsize=7)
par(mar=c(2.8,5.5,0.6,0.4),mgp=c(1.8,0.5,0),tcl=-0.3,lwd=.75,cex=1)
barplot(ex[6:1],horiz = T,log="x",border = NA,xlim=c(1,1e7),
        col=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.7,0.7,0.6,0.8,0.8))[6:1],
        names.arg = c("regular 2","ultra-clean 1","ultra-clean 2",
                      "regular 3","ultra-clean 3","ultra-clean 4")[6:1],las=1,
        xlab=expression(paste("potential contaminants [molecules / ", mu,"l]")) ,cex.lab=9/7)
dev.off()

#plot Venns - use web service for prettier version
pdf("data_for_Figure4_B.pdf",width=8.8/2.54,height=8.8/2.54,pointsize = 7)
anna3tab2venn(consTab[,c("A","B","C")],noCount = F)
dev.off()
pdf("data_for_Figure4_C.pdf",width=8.8/2.54,height=8.8/2.54,pointsize = 7)
anna3tab2venn(consTab[,c("D","E","F")],noCount = F)
dev.off()


#potential exogenous: not mapping to the human genome, not present on the columns, not present in the water sample
# but found in the plasma samples
shouldnots <- c(2:4,12:16,24:32,40:48,56:62) #where they should not be
shoulds <- c(5:11,17:23,33:39,49:55) #where they should be

clustBs <- clustB2[!rownames(clustB2) %in% c(clusAlg,clus8182,clus3378,clus3869,clus8030,clus6846),c(shouldnots,shoulds)]
clustBse <- clustBs[apply(clustBs[,1:33],1,function(x) length(which(x>0))<=3)&
                      apply(clustBs[,1:33],1,function(x) max(x)<=10)&
                      apply(clustBs[,34:61],1,function(x) length(which(x>3))>=8)&
                      rowMeans(clustBs[,34:61])>rowMeans(clustBs[,1:33]),]
clustBse <- clustBse[order(-rowSums(clustBse[,34:61])),]
potex <- clustB$Row.names[rownames(clustB) %in% rownames(clustBse)]
write.table(potex,"171018_potentialExogenousReads.txt",quote=F,row.names=F,col.names=F)

#sequences were extracted on the command line:
#for i in $(cat 171018_potentialExogenousReads.txt); do grep ^$i -w validUnmappedSequences.100.100.clustered.bak.clstr | grep -v %$ | cut -f 2 -d ">" | sed 's/[.][.][.]//g' | cut -f 1 -d " " | grep -f - -w -A 1 validUnmappedSequences.fa >> 171018_potentialExogenousReads.fa; done

#next, they were blasted and hits compared to the contaminants in Salter et al.
# clusters were identified on the command line
# for i in $(cat 171018_potentialExogenousReads.noHuman.noSalter.contignames.txt); do grep $i validUnmappedSequences.100.100.clustered.bak.clstr | cut -f 1 >> 171018_potentialExogenousReads.noHuman.noSalter.clusters; done
potc <- read.delim("171018_potentialExogenousReads.noHuman.noSalter.clusters",header=F,stringsAsFactors = F)

clust1e <- clustB[,-c(1,64,65)]
for(i in 1:ncol(clust1e)) clust1e[,i] <- clust1e[,i]/sum(clustT[,i])
write.table(clust1e[clustB$Row.names %in% potc$V1,]*1e6,"171031_potentialExogenousReads.noHuman.noSalter.clusters.cpm.tsv",
            sep="\t")
write.table(clustB[clustB$Row.names %in% potc$V1,"name"],"171031_potentialExogenousReads.noHuman.noSalter.clusters.cpm.tsv",
            sep="\t",append=T)
