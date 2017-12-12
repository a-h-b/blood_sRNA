library(scales)
library(gplots)

#read data from Excel and some re-formatting
ucpO <- read.delim("161013_qiaSummary_withDnase_extractUCTvsOld.txt",stringsAsFactors = F)
ucpO <- ucpO[order(ucpO$value,ucpO$plotname),]
ucpO <- data.frame(ucpO[,1:6],"space"=rep(0,nrow(ucpO)),ucpO[,7:9],"space2"=rep(0,nrow(ucpO)),stringsAsFactors = F)

pdf("Figure4A.pdf",width=8/2.54,height=5/2.54,pointsize=7)
par(mar=c(4.3,4.2,0.6,0.4),mgp=c(2.9,0.5,0),tcl=-0.3,lwd=.75,cex=1)
barplot2(t(ucpO[ucpO$value=="mean",c("regular.2.1","ultra.clean.1.1","ultra.clean.1.2","space","regular.2.2","ultra.clean.2.1","ultra.clean.2.2","space2")]),
         beside=T,border=NA,xlim=c(2.4,52),
         ylim=c(1,2E7),log="y",col=alpha(colors()[rep(c(645,116),each=4)],c(0.45,0.7,0.7,0.8,0.6,0.8,0.8,0.7)),xpd=F,
         names=ucpO$plotname[ucpO$value=="mean"],las=2,cex.names=9/7,ylab=expression(paste("molecules / ", mu,"l eluate")),cex.lab=9/7,
         axes=F,cex.axis=1,cex=1
         ,plot.ci=T,ci.l=ifelse(t(ucpO[ucpO$value=="mean",c("regular.2.1","ultra.clean.1.1","ultra.clean.1.2","space","regular.2.2","ultra.clean.2.1","ultra.clean.2.2","space2")])-
                                  t(ucpO[ucpO$value=="sd",c("regular.2.1","ultra.clean.1.1","ultra.clean.1.2","space","regular.2.2","ultra.clean.2.1","ultra.clean.2.2","space2")])>=1,
                                t(ucpO[ucpO$value=="mean",c("regular.2.1","ultra.clean.1.1","ultra.clean.1.2","space","regular.2.2","ultra.clean.2.1","ultra.clean.2.2","space2")])-
                                  t(ucpO[ucpO$value=="sd",c("regular.2.1","ultra.clean.1.1","ultra.clean.1.2","space","regular.2.2","ultra.clean.2.1","ultra.clean.2.2","space2")]),0.1),
         ci.u=t(ucpO[ucpO$value=="mean",c("regular.2.1","ultra.clean.1.1","ultra.clean.1.2","space","regular.2.2","ultra.clean.2.1","ultra.clean.2.2","space2")])+
           t(ucpO[ucpO$value=="sd",c("regular.2.1","ultra.clean.1.1","ultra.clean.1.2","space","regular.2.2","ultra.clean.2.1","ultra.clean.2.2","space2")]),
)
axis(2,las=2,yaxs="i",cex=1)
box(bty="l",lwd=1)
legend("top",legend=c("regular 2","ultra-clean 1&2","regular 3","ultra-clean 3&4"),border=NA,
       fil=alpha(colors()[rep(c(641,119),each=2)],c(0.45,0.7,0.6,0.8)),bty="n",
       title="eluate of column batch",ncol=2,x.intersp = 0.7)
dev.off()
