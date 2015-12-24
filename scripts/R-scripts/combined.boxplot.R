require("RColorBrewer")


#tools=c("RIsearch","IntaRNA","RNAplex","RNAcofold","pairfold","RNAup","RNAduplex","RNAhybrid","bifold","DuplexFold","ssearch","ractip","bistarna","AccessFold")

tools=c("RNAup","RNAplex","IntaRNA","RNAcofold","pairfold","DuplexFold","RNAduplex","RNAhybrid","ssearch","RIsearch","AccessFold","bifold","bistarna","ractip")


#Colors=rainbow(14)
Colors=c("#C3873B","#CB52CD","#769FBD","#6CD056","#583A2D","#C1467E","#C0C78C","#C6493A","#C7CC43","#4F6E3C","#50395F","#8672CA","#CD9CA7","#77D0B6")



bacterial.regional.results=read.csv("/home/suu13/projects/benchmark/bacteria/bacterial.results.table.regional.csv",header = F,sep="\t")
eukaryotic.regional.results=read.csv("/home/suu13/projects/benchmark/eukaryotes/results/eukaryotic.results.table.regional.csv",header = F,sep="\t")
archaeal.regional.results=read.csv("/home/suu13/projects/benchmark/archaea/results/archaeal.results.table.regional.csv",header = F,sep="\t")

extract_tool_result=function(tool,df){

  return(df[df$V1 == tool,])
  
}





results_for_boxplot=sapply(tools,extract_tool_result,df=rbind(bacterial.regional.results,eukaryotic.regional.results,archaeal.regional.results),USE.NAMES = T,simplify = F)

pdf("overall.results.pdf",width = 12.1,height = 11)
par(oma=c(10,4,2,2),mar=c(0,3,2,1),xpd=NA,cex.axis=1.1) #c(bottom, left, top, right)
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

#mtext("Overall performance of RNA interaction algorithms",outer = T,cex=1.8)



dev.off()