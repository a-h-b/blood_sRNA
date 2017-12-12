load("figure2_WS.Rdata")
library(gplots)

pdf("Figure2.pdf",width=12/2.54,height=10/2.54,pointsize=7)
layout(matrix(1:4,ncol=2,nrow=2,byrow=T))
par(mar=c(5,3.2,0.9,0.9),mgp=c(2,0.5,0),tcl=-0.3,lwd=.75,cex=1)
barplot2(t(f1val),beside=T,ylim=c(0,45),plot.ci=T,ci.l=t(f1val-f1sd),ci.u=t(f1val+f1sd),las=2,
         cex.names=9/7,col=colors()[c(76,516)],cex.lab=9/7,ylab="45 - Cp",axes=F,xpd=F,cex.axis=1,cex=1)
axis(2,las=2,yaxs="i",at=0:9*5,cex=1)
box(bty="l")
legend("topleft",legend=colnames(f1val),fil=colors()[c(76,516)],bty="n",
       title="input to cDNA synthesis")
barplot2(t(f2val),beside=T,ylim=c(0,45),plot.ci=T,ci.l=t(f2val-f2sd),ci.u=t(f2val+f2sd),las=2,
         cex.names=9/7,ylab="45 - Cp",col=colors()[c(430,116)],cex.lab=9/7,
         axes=F,cex.axis=1,cex=1)
axis(2,las=2,yaxs="i",at=0:9*5,cex=1)
box(bty="l")
legend("topleft",legend=gsub("mock ","mock-",colnames(f2val)),fil=colors()[c(430,116)],bty="n",
       title="input to cDNA synthesis")

par(mar=c(5,4.2,0.9,0.9),mgp=c(3,0.5,0),tcl=-0.3,lwd=.75,cex=1)
barplot2(t(fdnase[fdnase$value=="mean",c("mock","Dnase")]),beside=T,
         ylim=c(10,5E7),log="y",col=colors()[c(430,438)],xpd=F,
         names=fdnase$plotname[fdnase$value=="mean"],las=2,cex.names=9/7,ylab=expression(paste("molecules / ", mu,"l eluate")),cex.lab=9/7,
         axes=F,cex.axis=1,cex=1
         ,plot.ci=T,ci.l=ifelse(t(fdnase[fdnase$value=="mean",c("mock","Dnase")])-
                                  t(fdnase[fdnase$value=="sd",c("mock","Dnase")])>=1,t(fdnase[fdnase$value=="mean",c("mock","Dnase")])-
                                  t(fdnase[fdnase$value=="sd",c("mock","Dnase")]),0.1),
         ci.u=t(fdnase[fdnase$value=="mean",c("mock","Dnase")])+
           t(fdnase[fdnase$value=="sd",c("mock","Dnase")]),
)
axis(2,las=2,yaxs="i",cex=1)
box(bty="l")
legend("topleft",legend=c("mock-extract","+ DNase treatment"),fil=colors()[c(430,438)],bty="n",
       title="input to cDNA synthesis")

par(mar=c(5,3.2,0.9,0.9),mgp=c(2,0.5,0),tcl=-0.3,lwd=.75,cex=1)
barplot2(t(f3val),beside=T,ylim=c(0,100),plot.ci=T,ci.l=t(f3val-f3sd),ci.u=t(f3val+f3sd),las=2,xpd = F,
         cex.names=9/7,ylab="% sRNA remaining",col=c("white",colors()[129]),cex.lab=9/7,cex.axis=1,cex=1)
legend("topleft",legend=colnames(f3val),fil=c("white",colors()[129]),bty="n",
       title="pre-treatment of column")
box(bty="l")
dev.off()
