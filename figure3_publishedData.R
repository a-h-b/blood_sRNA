library(gplots)
library(RColorBrewer)
library(gtools)
source("heatmap.2m.R")

#read meta data and numbers of reads at different steps
metaX <- read.delim("allLibsMeta.txt",stringsAsFactors = F) #metadata
pub <- read.delim("public.readCountsPerSeqOIPerSample.unmapped.14.t.tsv",stringsAsFactors = F,row.names=1) #clustering result
tot <- read.delim("allReads.tsv",stringsAsFactors = F,row.names=1,header=F) #all reads

#merge read numbers
pubm <- merge(tot,pub,by.x=0,by.y=0,all.x=T)
for(i in 4:10) pubm[is.na(pubm[,i]),i] <- 0
#normalize reads to total number of reads
pubt <- as.matrix(pubm[,4:9])
rownames(pubt) <- pubm$Row.names
for(i in 1:6) pubt[,i] <- pubt[,i]/pubm$V3

#name missing metadata and add extraction kits
metaxt <- merge(pubt,metaX,by.x=0,by.y="Run_s",all=T)
metaxt$BioProject_s[grep("^SM",metaxt$Row.names)] <- "Wang et al."
metaxt$extraction <- ""
class(metaxt$extraction) <- "character"
metaxt$extraction[metaxt$BioProject_s %in% c("Wang et al.","PRJNA196121","PRJEB6683","PRJNA230622",
                                             "PRJNA213914","PRJNA268092")] <- "Q"
metaxt$extraction[metaxt$BioProject_s %in% c("PRJNA124945","PRJNA142489","PRJNA144887",
                                             "PRJNA140779","PRJNA127103","PRJNA112823","PRJNA116139",
                                             "PRJNA120351","PRJNA124019")] <- "T"
metaxt$extraction[metaxt$BioProject_s %in% c("PRJNA145935","PRJNA184379")] <- "V"
metaxt$extraction[metaxt$BioProject_s %in% c("PRJNA129787")] <- "P"

#per study data
plotMean0t <- aggregate(rowSums(metaxt[,2:7]),list(metaxt$BioProject_s),function(x) mean(x,na.rm=T))[,2]
plotSd0t <- aggregate(rowSums(metaxt[,2:7]),list(metaxt$BioProject_s),function(x) sd(x,na.rm=T))[,2]
plotNames0t <- aggregate(rowSums(metaxt[,2:7]),list(metaxt$BioProject_s),function(x) mean(x,na.rm=T))[,1]
plotNum0t <- aggregate(rowSums(metaxt[,2:7]),list(metaxt$BioProject_s),function(x) length(which(!is.na(x))))[,2]
plotEx0t <- as.character(aggregate(metaxt$extraction,list(metaxt$BioProject_s),function(x) unique(x))[,2])

#remove NAs
plotSd0t <- plotSd0t[!is.na(plotMean0t)]
plotNames0t <- plotNames0t[!is.na(plotMean0t)]
plotNum0t <- plotNum0t[!is.na(plotMean0t)]
plotEx0t <- plotEx0t[!is.na(plotMean0t)]
plotMean0t <- plotMean0t[!is.na(plotMean0t)]

#sort by mean
plotSd0t <- plotSd0t[order(plotMean0t,decreasing = T)]
plotNames0t <- plotNames0t[order(plotMean0t,decreasing = T)]
plotNum0t <- plotNum0t[order(plotMean0t,decreasing = T)]
plotEx0t <- plotEx0t[order(plotMean0t,decreasing = T)]
plotMean0t <- plotMean0t[order(plotMean0t,decreasing = T)]

#plot
pdf("Figure3.pdf",width=8/2.54,height=5/2.54,pointsize=7)
par(mar=c(6.3,4.2,2.6,0.4),mgp=c(3,0.5,0),tcl=-0.3,lwd=.75,cex=1)
a <- barplot2(1e6*plotMean0t,plot.ci = T,ci.l=1e6*ifelse(plotMean0t-plotSd0t>0&!is.na(plotMean0t-plotSd0t),plotMean0t-plotSd0t,1e-8),
              ci.u=1e6*(plotMean0t+plotSd0t),las=2,names.arg = plotNames0t,log="y",ylim=c(1,1e4),xpd=F,
              ylab="contaminant rpm",cex.lab=9/7,xlim=c(-0.8,22))
text(c(-0.6,a),1e5,labels = c("n =",plotNum0t),xpd=T)
text(c(-0.6,a),3e4,labels = c("E :",plotEx0t),xpd=T)
box(bty="l",lwd=1)
dev.off()


#add metadata for on type of material (from source, tissue and manually from articles)
plotmxt <- metaxt[metaxt$BioProject_s %in% plotNames0t[1:7],]
plotmxt$namenum <- sapply(plotmxt$BioProject_s,function(x) which(x %in% plotNames0t[1:7]))
plotmxt <- plotmxt[order(plotmxt$namenum,-rowSums(plotmxt[,2:7])),]
plotmxt$type <- unlist(sapply(gsub(";$","",gsub("^;","",tolower(paste(plotmxt$source_name_s,plotmxt$tissue_s,sep=";")))),
                              function(x) if(x=="") "" else unique(unlist(strsplit(x,split=";")))))
plotmxt$type[is.na(plotmxt$type)] <- ""
plotmxt$type[plotmxt$type=="na"] <- ""
plotmxt$type[plotmxt$BioProject_s=="PRJEB6683"] <- "plasma"
plotmxt$type[plotmxt$BioProject_s=="PRJNA213914"] <- "plasma"
plotmxt$type[grep("^SM",plotmxt$Row.names)] <- "plasma"

#prepare data for plotting
plotmxt2 <- as.matrix(plotmxt[,2:7])
rownames(plotmxt2) <- plotmxt$Row.names
rownames(plotmxt2) <- paste(plotmxt$BioProject_s,rownames(plotmxt2),sep="_")
plotmxt2 <- plotmxt2[rowSums(plotmxt2)>0,]
colnames(plotmxt2) <- c("sRNA 6","sRNA 2","sRNA 1","sRNA 5","sRNA 3","sRNA 4")
plotmxt2 <- plotmxt2[,order(colnames(plotmxt2))]


#elements for supplementary figure:

#heatmap with side colours
pdf("SupplementToFigure3.pdf",width=8/2.54,height=5/2.54,pointsize=6)
par(lwd=0.75)
hmret <- heatmap.2m(t(plotmxt2),scale="col",Rowv=NA,
                    Colv=as.dendrogram(hclust(as.dist(1-cor(t(plotmxt2),method="spearman")),method="ward.D2")),
                    col=colorRampPalette(brewer.pal(11,"RdYlBu"))(256)[256:1],density.info = "none",dendrogram="col",trace = "none",keyName = "contaminant",
                    ColSideColors = rbind(colorRampPalette(brewer.pal(9,"Greys"))(256)[cut(sapply(rownames(plotmxt2),
                                                                                                  function(x)sum(plotmxt[plotmxt$Row.names==gsub(".+_","",x),2:7])),256)],
                                          c("magenta","darkblue")[as.numeric(as.factor(sapply(rownames(plotmxt2),
                                                                                              function(x)plotmxt$extraction[plotmxt$Row.names==gsub(".+_","",x)])))],
                                          brewer.pal(7,"Set2")[as.numeric(as.factor(sapply(rownames(plotmxt2),
                                                                                           function(x)plotmxt$type[plotmxt$Row.names==gsub(".+_","",x)])))],
                                          brewer.pal(7,"Pastel2")[as.numeric(as.factor(gsub("_.+","",rownames(plotmxt2))))]
                                          
                    ),
                    margins = c(7,5),ColSideFac = 1.6,cexCol = 0.6,cexRow = 0.8)
dev.off()

#colour bars:
pdf("SupplementToFigure3_colorBar.pdf",width=0.55/2.54,height=2.4/2.54,pointsize=4)
par(mar=c(0.25,1.35,0.25,0.25),tcl=-0.3,mgp=c(1.5,0.45,0))
image(t(cbind(1:length(hmret$col),1:length(hmret$col))),
      col=hmret$col,axes=F)
axis(2,at=seq(from=1,to=256,length.out = 5)/256,las=1,
     labels=round(hmret$breaks[round(c(seq(from=round(256*1/18),to=round(256*18/18),length.out = 5)))]),cex.axis=0.8)
box()
dev.off()

pdf("SupplementToFigure3_contaminantBar.pdf",width=2.2/2.54,height=0.6/2.54,pointsize=4)
par(mar=c(1.3,0.75,0.5,0.75),tcl=-0.3,mgp=c(1.5,0.25,0))
image(cbind(1:256,1:256),col=colorRampPalette(brewer.pal(9,"Greys"))(256),axes=F)
axis(1,at=c(2,128,255)/256,
     labels=round(1e6*as.numeric(gsub("\\(","",gsub(",.+","",levels(cut(sapply(rownames(plotmxt2),
                                                                               function(x)sum(plotmxt[plotmxt$Row.names==gsub(".+_","",x),2:7])),256))[c(2,128,255)])))),cex.axis=0.8)
box()
dev.off()

#legends:
pdf("SupplementToFigure3_legend.pdf",width=10/2.54,height=15/2.54,pointsize=6)
plot(1,type="n",axes=F,ann=F)
legend("topright",legend=levels(as.factor(sapply(rownames(plotmxt2),
                                                 function(x)plotmxt$extraction[plotmxt$Row.names==gsub(".+_","",x)]))),
       fill = c("magenta","darkblue"),bty="n",border = F)
legend("bottomright",legend=gsub("mouse primary ","",levels(as.factor(sapply(rownames(plotmxt2),
                                                                             function(x)plotmxt$type[plotmxt$Row.names==gsub(".+_","",x)])))),
       fill = brewer.pal(7,"Set2")[1:4],bty="n",border = F)
legend("topleft",legend=levels(as.factor(gsub("_.+","",rownames(plotmxt2)))),
       fill =  brewer.pal(7,"Pastel2"),bty="n",border = F)
dev.off()



