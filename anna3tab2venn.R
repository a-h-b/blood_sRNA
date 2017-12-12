library(RColorBrewer)

anna3tab2venn <- function(data,noCount="unknown",col=brewer.pal(3,"Set1")){
  require(VennDiagram)
  require(scales)
  plot(1,type="n",ann=F,axes=F)
  draw.triple.venn(area1 =length(which(data[,1] != noCount)),
                   area2 =length(which(data[,2] != noCount)),
                   area3 =length(which(data[,3] != noCount)),
                   n12 =length(which(data[,1] != noCount & data[,2] != noCount)), 
                   n13=length(which(data[,1] != noCount & data[,3] != noCount)), 
                   n23=length(which(data[,2] != noCount & data[,3] != noCount)), 
                   n123=length(which(data[,1] != noCount & data[,2] != noCount& data[,3] != noCount)), 
                   category=colnames(data),col=col,fil=col,alpha=rep(0.3,3),cex=0.8,fontfamily=rep("sans",7),
                   cat.fontfamily=rep("sans",3)
  )
}
