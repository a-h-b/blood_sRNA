library(scales)
library(gplots)

#read data
ucpO <- read.delim("matchcountSeqOIPerSample4plot.txt",comment.char = "#",stringsAsFactors = F)
ucpOp <- ucpO
ucpOp[ucpOp==0] <- 0.001 #pseudocount for plotting
colnames(ucpO) <- gsub("regular.2.1","regular.2",
                       gsub("ultra.clean.1.1","ultra.clean.1",
                            gsub("ultra.clean.1.2","ultra.clean.2",
                                 gsub("regular.2.2","regular.3",
                                      gsub("ultra.clean.2.1","ultra.clean.3",
                                           gsub("ultra.clean.2.2","ultra.clean.4",
                                                colnames(ucpO)))))))


pdf("Supplementary_Figure3.pdf",width=13/2.54,height=13/2.54,pointsize=7)
layout(matrix(1:8,ncol=4,nrow=2,byrow=T))
par(mar=c(7,4.8,0.6,0.4),mgp=c(3.5,0.4,0),tcl=-0.3,lwd=.75,cex=1)
for(i in 1:6){
  ylaba <- paste0("sRNA",i," [molecules /")
  ylabb <- "l eluate]"
  barplot2(colSums(as.matrix(ucpO[10:13,-1][,6*(i-1)+1:6])),
           #beside=F,
           border=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.7,0.7,0.5,0.8,0.8)),#xlim=c(2.4,52),
           col="white",lwd=1.5,
           ylim=c(0,1.1*max(colSums(as.matrix(ucpO[10:13,-1][,6*(i-1)+1:6])))),#log="y",
           xpd=F,
           names=gsub("regular.","regular ",gsub("ultra.clean.","ultra-clean ",gsub("sRNA1.","",colnames(ucpO)[2:7]))),las=2,cex.names=1,
           ylab=bquote(.(ylaba) ~ mu*.(ylabb)),
           cex.lab=9/7,
           axes=F,cex.axis=1,cex=1
           ,plot.ci=F)
  barplot2(colSums(as.matrix(ucpO[10:13,-1][,6*(i-1)+1:6])),add=T,
           #beside=F,
           border=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.7,0.7,0.5,0.8,0.8)),#xlim=c(2.4,52),
           col=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.7,0.7,0.5,0.8,0.8)),density=5,angle=135,
           ylim=c(0,1.1*max(colSums(as.matrix(ucpO[10:13,-1][,6*(i-1)+1:6])))),#log="y",
           xpd=F,
           names=rep("",6),#las=2,cex.names=9/7,
           ylab="",cex.lab=9/7,
           axes=F,cex.axis=1,cex=1
           ,plot.ci=F)
  barplot2(colSums(as.matrix(ucpO[10:12,-1][,6*(i-1)+1:6])),add = T,
           #beside=F,
           border=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.7,0.7,0.5,0.8,0.8)),#xlim=c(2.4,52),
           col="white",
           ylim=c(0,1.1*max(colSums(as.matrix(ucpO[10:13,-1][,6*(i-1)+1:6])))),#log="y",
           xpd=F,
           names=rep("",6),#las=2,cex.names=9/7,
           ylab="",cex.lab=9/7,
           axes=F,cex.axis=1,cex=1
           ,plot.ci=F  )
  barplot2(colSums(as.matrix(ucpO[10:12,-1][,6*(i-1)+1:6])),add = T,
           #beside=F,
           border=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.7,0.7,0.5,0.8,0.8)),#xlim=c(2.4,52),
           col=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.7,0.7,0.5,0.8,0.8)),density=45,
           ylim=c(0,1.1*max(colSums(as.matrix(ucpO[10:13,-1][,6*(i-1)+1:6])))),#log="y",
           xpd=F,
           names=rep("",6),#las=2,cex.names=9/7,
           ylab="",cex.lab=9/7,
           axes=F,cex.axis=1,cex=1
           ,plot.ci=F  )
  barplot2(colSums(as.matrix(ucpO[10:11,-1][,6*(i-1)+1:6])),add = T,
           #beside=F,
           border=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.7,0.7,0.5,0.8,0.8)),#xlim=c(2.4,52),
           col="white",
           ylim=c(0,1.1*max(colSums(as.matrix(ucpO[10:13,-1][,6*(i-1)+1:6])))),#log="y",
           xpd=F,
           names=rep("",6),#las=2,cex.names=9/7,
           ylab="",cex.lab=9/7,
           axes=F,cex.axis=1,cex=1
           ,plot.ci=F  )
  barplot2(colSums(as.matrix(ucpO[10:11,-1][,6*(i-1)+1:6])),add = T,
           #beside=F,
           border=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.7,0.7,0.5,0.8,0.8)),#xlim=c(2.4,52),
           col=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.7,0.7,0.5,0.8,0.8)),density=60,
           ylim=c(0,1.1*max(colSums(as.matrix(ucpO[10:13,-1][,6*(i-1)+1:6])))),#log="y",
           xpd=F,
           names=rep("",6),#las=2,cex.names=9/7,
           ylab="",cex.lab=9/7,
           axes=F,cex.axis=1,cex=1
           ,plot.ci=F  )
  barplot2(as.numeric(ucpO[10,-1][6*(i-1)+1:6]),add = T,
           #beside=F,
           border=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.7,0.7,0.5,0.8,0.8)),#xlim=c(2.4,52),
           col="white",
           ylim=c(0,1.1*max(colSums(as.matrix(ucpO[10:13,-1][,6*(i-1)+1:6])))),#log="y",
           xpd=F,
           names=rep("",6),#las=2,cex.names=9/7,
           ylab="",cex.lab=9/7,
           axes=F,cex.axis=1,cex=1
           ,plot.ci=F  )
  barplot2(as.numeric(ucpO[10,-1][6*(i-1)+1:6]),add = T,
           #beside=F,
           border=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.7,0.7,0.5,0.8,0.8)),#xlim=c(2.4,52),
           col=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.7,0.7,0.5,0.8,0.8)),
           ylim=c(0,1.1*max(colSums(as.matrix(ucpO[10:13,-1][,6*(i-1)+1:6])))),#log="y",
           xpd=F,
           names=rep("",6),#las=2,cex.names=9/7,
           ylab="",cex.lab=9/7,
           axes=F,cex.axis=1,cex=1
           ,plot.ci=F  )
  pts <- pretty(c(0,max(colSums(as.matrix(ucpO[10:14,-1][,6*(i-1)+1:6])))))
  axis(2, at = pts, labels = format(pts,scientific = T,digits = 2),las=2)
  box(bty="l",lwd=1)
}
par(mar=c(0,0,0,0))
plot(1,type="n",ann=F,axes=F)
legend("center",legend=c("1","2",">2"),border="grey30",
       fill ="grey30",density = c(40,60,105),bty="n",cex=1,
       title="mismatches to\nhuman genome",x.intersp = 0.7,ncol=2)
plot(1,type="n",ann=F,axes=F)
legend("center",legend=c("regular 2","ultra-clean 1&2","regular 3","ultra-clean 3&4"),border=NA,
       fill =alpha(colors()[rep(c(645,116),each=2)],c(0.45,0.75,0.5,0.8)),
       bty="n",cex=1,
       title="eluate of column batch",x.intersp = 0.7,ncol=1)
dev.off()

