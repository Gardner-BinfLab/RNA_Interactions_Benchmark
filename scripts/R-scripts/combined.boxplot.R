require("RColorBrewer")

require("fields")
require("gplots")

require("plotrix")


#tools=c("RIsearch","IntaRNA","RNAplex","RNAcofold","pairfold","RNAup","RNAduplex","RNAhybrid","bifold","DuplexFold","ssearch","ractip","bistarna","AccessFold")

tools=c("RNAup","IntaRNA","RNAplex","RIsearch","RNAcofold","pairfold","DuplexFold","RNAduplex","RNAhybrid","ssearch","AccessFold","bifold","bistaRNA","ractIP")


#Colors=rainbow(14)
Colors=c("#C3873B","#CB52CD","#769FBD","#6CD056","#583A2D","#C1467E","#C0C78C","#C6493A","#C7CC43","#4F6E3C","#50395F","#8672CA","#CD9CA7","#77D0B6")



bacterial.regional.results=read.csv("/home/suu13/projects/benchmark/bacteria/bacterial.results.table.regional.csv",header = F,sep="\t")
eukaryotic.regional.results=read.csv("/home/suu13/projects/benchmark/eukaryotes/results/eukaryotic.results.table.regional.csv",header = F,sep="\t")
archaeal.regional.results=read.csv("/home/suu13/projects/benchmark/archaea/results/archaeal.results.table.regional.csv",header = F,sep="\t")

#bacterial.regional.results=read.csv("bacterial.results.table.regional.csv",header = F,sep="\t")
#eukaryotic.regional.results=read.csv("eukaryotic.results.table.regional.csv",header = F,sep="\t")
#archaeal.regional.results=read.csv("archaeal.results.table.regional.csv",header = F,sep="\t")


extract_tool_result=function(tool,df){

  return(df[tolower(df$V1) == tolower(tool),])
  
}





results_for_boxplot=sapply(tools,extract_tool_result,df=rbind(bacterial.regional.results,eukaryotic.regional.results,archaeal.regional.results),USE.NAMES = T,simplify = F)

pdf("overall.results.boxplot.pdf",width = 12.1,height = 11)
par(oma=c(10,4,2,2),mar=c(1,3,2,1),xpd=NA,cex.axis=1.1) #c(bottom, left, top, right)
layout(matrix(c(1,1,2,2,3,3),nrow=3,byrow = T))
boxplot(do.call(cbind,lapply(results_for_boxplot,function(df) return(df$V5))),names=F,col=Colors,ylim=c(-0.05,1.05)) #TPR
mtext("A",side = 3,adj=0,cex=1.8,font=2)
mtext("TPR (Sensitivity)",side=2,line=3,cex=1.8)
boxplot(do.call(cbind,lapply(results_for_boxplot,function(df) return(df$V6))),names=F,col=Colors,ylim=c(-0.05,1.05)) #PPV
mtext("B",side = 3,adj=0,cex=1.8,font=2)
mtext("PPV (Precision)",side=2,line=3,cex=1.8)
boxplot(do.call(cbind,lapply(results_for_boxplot,function(df) return(df$V7))),names=F,col=Colors,ylim=c(-0.05,1.05)) #MCC
mtext("C",side = 3,adj=0,cex=1.8,font=2)
mtext("MCC",side=2,line=3,cex=1.8)
axis(1,at=1:14,cex.axis=1.8,las=2,labels=F)
text(1:14,-0.3,labels=tools,cex=2,srt=45)

dev.off()



cor.matrix.tools=cor(do.call(cbind,lapply(results_for_boxplot,function(df) return(df$V7))))
cor.matrix.tools.reorder=cor.matrix.tools[hclust(dist(cor.matrix.tools))$order,hclust(dist(t(cor.matrix.tools)))$order]

pdf("overall.results.correlation.pdf",width = 12.1*1.5,height = 11*1.5)
par(mar=c(10,10,2,7),xpd=NA,cex=1.8) #c(bottom, left, top, right)
image(cor.matrix.tools.reorder,col=two.colors(n=12,"cadetblue","firebrick","grey"),axes=F)

axis(2,at=seq(0,1,by=1/13),labels = F,cex=1.8)
text(-0.2,seq(0,1,by=1/13),labels=rownames(cor.matrix.tools.reorder),cex=1.5,srt=0,xpd=T)

axis(1,at=seq(0,1,by=1/13),labels = colnames(cor.matrix.tools.reorder),las=2,cex.axis=1.5)
#text(seq(0,1,by=1/13),-0.15,labels=colnames(cor.matrix.tools.reorder),cex=1.5,srt=90,xpd=T)

image.plot(cor.matrix.tools.reorder,col=two.colors(n=12,"cadetblue","firebrick","grey"),axes=F,legend.only=T,axis.args = list(cex.axis = 1.8,las=1))


dev.off()








pdf("overall.results.pheatmap.pdf",onefile=FALSE,width = 10)
pheatmap(tools_result_matrix,cluster_rows=FALSE,width = 3,fontsize_col=12)


dev.off()



pdf("overall.results.heatmap.pdf",width = 12.1,height = 11*1.5)

tools_result_matrix=do.call(cbind,lapply(results_for_boxplot,function(df) return(df$V7)))
heatmap.2(tools_result_matrix,Rowv = F,margins = c(10,10),col=colorRampPalette(rev(brewer.pal(n = 7, name =                                                                                                "RdYlBu")))(100),
          cexCol = 1.8, lhei = c(1.5, 1),trace="none",key=F,keysize=1,key.title ="",density.info="density",cexRow = 0.8,labRow=as.character(lapply(results_for_boxplot,function(df) return(df$V8))$RIsearch))

par(fig=c(0.2,0.95,0.8,1),xpd=T,mar=c(0,0,0,0)) #c(bottom, left, top, right) 
mtext("MCC",side=4,xpd=T,outer=F,line=-5,cex = 1.5)
image.plot(tools_result_matrix,legend.only = T,col=two.colors(n=100,"lightsteelblue4","firebrick","grey"),horizontal = F)
dev.off()



#more stats of predictions
aggregate(do.call(rbind,results_for_boxplot),by=list(do.call(rbind,results_for_boxplot)$V8),mean)

