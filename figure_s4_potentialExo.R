library(gplots)
library(scales)

exR <- read.delim("exoRNAs_values.txt",stringsAsFactors = F)
xcoord <- c(1,3:7)

pdf("SFigure4.pdf",width=29.5/2.54,height=6.1/2.54,pointsize=7)
par(mar=c(4.1,3.3,0.6,0.4),mgp=c(3,0.5,0),tcl=-0.3,lwd=0.75,cex=1)
plot(1,type="n",ann=F,axes=F,xlim=c(0,21*9),ylim=c(0,50),xaxs="i")
mtext("potential exogenous reads [cpm]",2,2.1,cex=9/7)
for(i in 1:nrow(exR)){
  lines(xcoord[2:6]+9*(i-1),
        exR[i,c("regular.2.45","regular.2.100","regular.2.223","regular.2.500","regular.2.1115")],
        col="darkgoldenrod3",lwd=1)
  lines(rep(xcoord[1]+9*(i-1),2),
        c(exR[i,"mean.neg"],sum(exR[i,c("mean.neg","sd.meg")])),
        col="grey",lwd=1)
  lines(rep(xcoord[3]+9*(i-1),2),
        c(exR[i,"regular.2.100"],sum(exR[i,c("regular.2.100","regular.2.100.sd")])),
        col="grey",lwd=1)
  points(xcoord+9*(i-1),
         exR[i,c("mean.neg","regular.2.45","regular.2.100","regular.2.223","regular.2.500","regular.2.1115")],
         col=c("lightskyblue3",rep("darkgoldenrod3",5)),pch=16)
}
for(i in 1:nrow(exR)){
  lines(xcoord[2:6]+9*(i-1),
        exR[i,c("ucp.1.45","ucp.1.100","ucp.1.223","ucp.1.500","ucp.1.1115")],
        col="darkgoldenrod1",lwd=1)
  lines(rep(xcoord[3]+9*(i-1),2),
        c(exR[i,"ucp.1.100"],sum(exR[i,c("ucp.1.100","ucp.1.100.sd")])),
        col="grey",lwd=1)
  points(xcoord[-1]+9*(i-1),
         exR[i,c("ucp.1.45","ucp.1.100","ucp.1.223","ucp.1.500","ucp.1.1115")],
         col=rep("darkgoldenrod1",5),pch=16)
}
for(i in 1:nrow(exR)){
  lines(xcoord[2:6]+9*(i-1),
        exR[i,c("regular.3.45","regular.3.100","regular.3.223","regular.3.500","regular.3.1115")],
        col="gold3",lwd=1)
  lines(rep(xcoord[3]+9*(i-1),2),
        c(exR[i,"regular.3.100"],sum(exR[i,c("regular.3.100","regular.3.100.sd")])),
        col="grey",lwd=1)
  points(xcoord[-1]+9*(i-1),
         exR[i,c("regular.3.45","regular.3.100","regular.3.223","regular.3.500","regular.3.1115")],
         col=rep("gold3",5),pch=16)
}
for(i in 1:nrow(exR)){
  lines(xcoord[2:6]+9*(i-1),
        exR[i,c("ucp.3.45","ucp.3.100","ucp.3.223","ucp.3.500","ucp.3.1115")],
        col="gold1",lwd=1)
  lines(rep(xcoord[3]+9*(i-1),2),
        c(exR[i,"ucp.3.100"],sum(exR[i,c("ucp.3.100","ucp.3.100.sd")])),
        col="grey",lwd=1)
  points(xcoord[-1]+9*(i-1),
         exR[i,c("ucp.3.45","ucp.3.100","ucp.3.223","ucp.3.500","ucp.3.1115")],
         col=rep("gold1",5),pch=16)
}
axis(1,at=c(0:20*9+4.1),labels=exR$X,mgp=c(3,2.8,0),tick=F,cex.axis=9/7,las=1)
axis(1,at=c(1:(9*21))[-c(rep(0:20*9,each=3)+c(2,8,9))],labels=rep(c(0,45,100,223,500,1115),times=21),las=2)
axis(2,las=2)
box(bty="l",lwd=1)
legend("top",legend="all controls",border=NA,col="lightskyblue3",pch=16,bty="n",title="")
legend("topright",legend=c("regular 2","ultra-clean 1","regular 3","ultra-clean 3"),border=NA,
       col=c("darkgoldenrod3","darkgoldenrod1","gold3","gold1"),bty="n",lty=1,lwd=1,pch=16,
       title="plasma extracted with column batch",ncol=5,x.intersp = 0.9)
dev.off()
