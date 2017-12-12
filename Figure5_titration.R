library(gplots)
library(scales)

ucpT <- read.delim("titration4plot.wostop.txt",stringsAsFactors = F)
ucpT$xcoord <- 7*(as.numeric(as.factor(ucpT$plotname))-1)+as.numeric(as.factor(ucpT$value))

pdf("Figure5A_titration.pdf",width=10.7/2.54,height=6.1/2.54,pointsize=7)
par(mar=c(4.1,4.2,0.6,0.4),mgp=c(3,0.5,0),tcl=-0.3,lwd=0.75,cex=1)
plot(1,type="n",ann=F,axes=F,xlim=c(1,max(ucpT$xcoord)),ylim=c(10,2e4),log="y")
mtext("contaminant reads [cpm]",2,3,cex=9/7)
for(i in 1:length(unique(ucpT$plotname))){ 
  lines(ucpT$xcoord[ucpT$plotname==unique(ucpT$plotname)[i]],
        ifelse(ucpT$regular.2.1[ucpT$plotname==unique(ucpT$plotname)[i]]>=1,
               ucpT$regular.2.1[ucpT$plotname==unique(ucpT$plotname)[i]],
               0.1),
        col=alpha(colors()[rep(c(645,116),each=2)],c(0.45,0.75,0.5,0.8))[1],lwd=1)
  points(ucpT$xcoord[ucpT$plotname==unique(ucpT$plotname)[i]],
         ifelse(ucpT$regular.2.1[ucpT$plotname==unique(ucpT$plotname)[i]]>=1,
                ucpT$regular.2.1[ucpT$plotname==unique(ucpT$plotname)[i]],
                0.1),
         col=alpha(colors()[rep(c(645,116),each=2)],c(0.45,0.75,0.5,0.8))[1],pch=16)
}
for(i in 1:length(unique(ucpT$plotname))){
  lines(ucpT$xcoord[ucpT$plotname==unique(ucpT$plotname)[i]],
        ifelse(ucpT$ultra.clean.1.1[ucpT$plotname==unique(ucpT$plotname)[i]]>=1,
               ucpT$ultra.clean.1.1[ucpT$plotname==unique(ucpT$plotname)[i]],
               0.1),
        col=alpha(colors()[rep(c(645,116),each=2)],c(0.45,0.75,0.5,0.8))[2],lwd=1)
  points(ucpT$xcoord[ucpT$plotname==unique(ucpT$plotname)[i]],
         ifelse(ucpT$ultra.clean.1.1[ucpT$plotname==unique(ucpT$plotname)[i]]>=1,
                ucpT$ultra.clean.1.1[ucpT$plotname==unique(ucpT$plotname)[i]],
                0.1),
         col=alpha(colors()[rep(c(645,116),each=2)],c(0.45,0.75,0.5,0.8))[2],pch=16)
}
for(i in 1:length(unique(ucpT$plotname))){
  lines(ucpT$xcoord[ucpT$plotname==unique(ucpT$plotname)[i]&!is.na(ucpT$regular.2.2)],
        ifelse(ucpT$regular.2.2[ucpT$plotname==unique(ucpT$plotname)[i]&!is.na(ucpT$regular.2.2)]>=1,
               ucpT$regular.2.2[ucpT$plotname==unique(ucpT$plotname)[i]&!is.na(ucpT$regular.2.2)],
               0.1),
        col=alpha(colors()[rep(c(645,116),each=2)],c(0.45,0.75,0.5,0.8))[3],lwd=1)
  points(ucpT$xcoord[ucpT$plotname==unique(ucpT$plotname)[i]&!is.na(ucpT$regular.2.2)],
         ifelse(ucpT$regular.2.2[ucpT$plotname==unique(ucpT$plotname)[i]&!is.na(ucpT$regular.2.2)]>=1,
                ucpT$regular.2.2[ucpT$plotname==unique(ucpT$plotname)[i]&!is.na(ucpT$regular.2.2)],
                0.1),
         col=alpha(colors()[rep(c(645,116),each=2)],c(0.45,0.75,0.5,0.8))[3],pch=16)
}
for(i in 1:length(unique(ucpT$plotname))){
  lines(ucpT$xcoord[ucpT$plotname==unique(ucpT$plotname)[i]&!is.na(ucpT$ultra.clean.2.1)],
        ifelse(ucpT$ultra.clean.2.1[ucpT$plotname==unique(ucpT$plotname)[i]&!is.na(ucpT$ultra.clean.2.1)]>=1,
               ucpT$ultra.clean.2.1[ucpT$plotname==unique(ucpT$plotname)[i]&!is.na(ucpT$ultra.clean.2.1)],
               0.1),
        col=alpha(colors()[rep(c(645,116),each=2)],c(0.45,0.75,0.5,0.8))[4],lwd=1)
  points(ucpT$xcoord[ucpT$plotname==unique(ucpT$plotname)[i]&!is.na(ucpT$ultra.clean.2.1)],
         ifelse(ucpT$ultra.clean.2.1[ucpT$plotname==unique(ucpT$plotname)[i]&!is.na(ucpT$ultra.clean.2.1)]>=1,
                ucpT$ultra.clean.2.1[ucpT$plotname==unique(ucpT$plotname)[i]&!is.na(ucpT$ultra.clean.2.1)],
                0.1),
         col=alpha(colors()[rep(c(645,116),each=2)],c(0.45,0.75,0.5,0.8))[4],pch=16)
}
abline(h=1e2,lwd=0.7,lty=3)
axis(1,at=c(1:41)[-c(7,14,21,28,35)],labels=rep(c(0,45,100,223,500,1115),times=6),las=2)
axis(1,at=c(0:5*7+3.6),labels=unique(ucpT$plotname),mgp=c(3,2.8,0),tick=F,cex.axis=9/7,las=1)
axis(2,las=2)
box(bty="l",lwd=1)
legend("top",legend=c("regular 2","ultra-clean 1","regular 3","ultra-clean 3"),border=NA,
       col=alpha(colors()[rep(c(645,116),each=2)],c(0.45,0.75,0.5,0.8)),bty="n",lty=1,lwd=1,pch=16,
       title="plasma extracted with column batch",ncol=4,x.intersp = 0.9)
dev.off()

ucp100 <- read.delim("readCountsPerSeqOIPerSample.unmapped.14.100forplot.txt",stringsAsFactors = F)
ucp100 <- ucp100[order(ucp100$value,ucp100$plotname),]
ucp100 <- data.frame(ucp100[,1:5],"space"=rep(0,nrow(ucp100)),ucp100[,6:7],"space2"=rep(0,nrow(ucp100)),stringsAsFactors = F)

pdf("Figure5B_100ulbatches.pdf",width=6/2.54,height=4.3/2.54,pointsize=7)
par(mar=c(4.3,4.2,0.6,0.4),mgp=c(3,0.5,0),tcl=-0.3,lwd=.75,cex=1)
barplot2(t(ucp100[ucp100$value=="mean",c("regular.2.1","ultra.clean.1.1","space","regular.2.2","ultra.clean.2.1","space2")]),
         beside=T,border=NA,xlim=c(0.4,40),
         ylim=c(1,2E5),log="y",col=alpha(colors()[rep(c(645,116),each=3)],c(0.45,0.75,0.8,0.5,0.8,0.7)),xpd=F,
         names=ucp100$plotname[ucp100$value=="mean"],las=2,cex.names=9/7,ylab="contaminants [cpm]",cex.lab=9/7,
         axes=F,cex.axis=1,cex=1
         ,plot.ci=T,ci.l=ifelse(t(ucp100[ucp100$value=="mean",c("regular.2.1","ultra.clean.1.1","space","regular.2.2","ultra.clean.2.1","space2")])-
                                  t(ucp100[ucp100$value=="sd",c("regular.2.1","ultra.clean.1.1","space","regular.2.2","ultra.clean.2.1","space2")])>=1,
                                t(ucp100[ucp100$value=="mean",c("regular.2.1","ultra.clean.1.1","space","regular.2.2","ultra.clean.2.1","space2")])-
                                  t(ucp100[ucp100$value=="sd",c("regular.2.1","ultra.clean.1.1","space","regular.2.2","ultra.clean.2.1","space2")]),0.1),
         ci.u=t(ucp100[ucp100$value=="mean",c("regular.2.1","ultra.clean.1.1","space","regular.2.2","ultra.clean.2.1","space2")])+
           t(ucp100[ucp100$value=="sd",c("regular.2.1","ultra.clean.1.1","space","regular.2.2","ultra.clean.2.1","space2")]),
)
axis(2,las=2,yaxs="i",cex=1)
box(bty="l",lwd=1)
legend("top",legend=c("regular 2","ultra-clean 1","regular 3","ultra-clean 3"),border=NA,
       fil=alpha(colors()[rep(c(645,116),each=2)],c(0.45,0.75,0.5,0.8)),bty="n",
       title="plasma extracted with column batch",ncol=2,x.intersp = 0.7)
dev.off()
