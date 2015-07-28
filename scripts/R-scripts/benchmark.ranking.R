

target_per_sRNA=function(number,tool.result.df){
  
  if(number== 0) return(0);
  total=0
  for (i in levels(tool.result.df$Qname) ) {
    total=total+as.numeric(table(tool.result.df[tool.result.df$Qname == i,]$TP[1:number])['TRUE'])
  }
  return(total)
}


mark_true_interactions=function(sRNA_locus_tag,mRNA_locus_tag){
  sRNA_info=bacterial.targets[bacterial.targets$Locus_Tag_.from.NCBI. == sRNA_locus_tag,]
  result=grepl(mRNA_locus_tag,paste0(sRNA_info$Target_Locus_Tag_3,sRNA_info$Target_Locus_Tag_2,sRNA_info$Target_Locus_Tag_1),fixed = T)
  return(result)
  
}
  
sRNA_TP_rank_find=function(sRNA_locus_tag,tool_DF){
  
  return(which(tool_DF[tool_DF$Qname ==sRNA_locus_tag,]$TP == TRUE))
}


bacterial.targets=read.table("/home/suu13/projects/benchmark/bacteria/bacterial.targets.tsv",sep="\t",header=T)

sRNA_locus_tags=levels(bacterial.targets$Locus_Tag_.from.NCBI.)

RIsearch.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/RIsearch/RIsearch.result.txt")
IntaRNA.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/IntaRNA/IntaRNA.result.csv")
RNAplex.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/RNAplex/RNAplex.result.csv")
RNAcofold.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/RNAcofold/RNAcofold.result.csv")
Pairfold.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/pairfold/pairfold.result.txt")
RNAup.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/RNAup/Native/RNAup.result.csv")
RNAduplex.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/RNAduplex/RNAduplex.result.csv")
RNAhybrid.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/RNAhybrid/RNAhybrid.result.csv")



names(RIsearch.results)=c("Qname","Qbeg","Qend", "Tname", "Tbeg", "Tend", "score", "energy")
names(IntaRNA.results)=c("Tname","Qname","energy")
names(RNAplex.results)=c("Tname","Qname","energy")
names(RNAcofold.results)=c("Qname","Tname","energy")
names(Pairfold.results)=c("Qname","Tname","energy")
names(RNAup.results)=c("Tname","Qname","energy")
names(RNAduplex.results)=c("Qname","Tname","energy")
names(RNAhybrid.results)=c("Tname","Qname","energy")

RIsearch.results=cbind(RIsearch.results,TP=as.factor(mapply(mark_true_interactions,as.character(RIsearch.results$Qname),as.character(RIsearch.results$Tname),USE.NAMES = F)))
IntaRNA.results=cbind(IntaRNA.results,TP=as.factor(mapply(mark_true_interactions,as.character(IntaRNA.results$Qname),as.character(IntaRNA.results$Tname),USE.NAMES = F)))
RNAplex.results=cbind(RNAplex.results,TP=as.factor(mapply(mark_true_interactions,as.character(RNAplex.results$Qname),as.character(RNAplex.results$Tname),USE.NAMES = F)))
RNAcofold.results=cbind(RNAcofold.results,TP=as.factor(mapply(mark_true_interactions,as.character(RNAcofold.results$Qname),as.character(RNAcofold.results$Tname),USE.NAMES = F)))
Pairfold.results=cbind(Pairfold.results,TP=as.factor(mapply(mark_true_interactions,as.character(Pairfold.results$Qname),as.character(Pairfold.results$Tname),USE.NAMES = F)))
RNAup.results=cbind(RNAup.results,TP=as.factor(mapply(mark_true_interactions,as.character(RNAup.results$Qname),as.character(RNAup.results$Tname),USE.NAMES = F)))
RNAduplex.results=cbind(RNAduplex.results,TP=as.factor(mapply(mark_true_interactions,as.character(RNAduplex.results$Qname),as.character(RNAduplex.results$Tname),USE.NAMES = F)))
RNAhybrid.results=cbind(RNAhybrid.results,TP=as.factor(mapply(mark_true_interactions,as.character(RNAhybrid.results$Qname),as.character(RNAhybrid.results$Tname),USE.NAMES = F)))



RIsearch.results=RIsearch.results[order(RIsearch.results$energy,decreasing = F),]
IntaRNA.results=IntaRNA.results[order(IntaRNA.results$energy,decreasing = F),]
RNAplex.results=RNAplex.results[order(RNAplex.results$energy,decreasing = F),]
RNAcofold.results=RNAcofold.results[order(RNAcofold.results$energy,decreasing = F),]
Pairfold.results=Pairfold.results[order(Pairfold.results$energy,decreasing = F),]
RNAup.results=RNAup.results[order(RNAup.results$energy,decreasing = F),]
RNAduplex.results=RNAduplex.results[order(RNAduplex.results$energy,decreasing = F),]
RNAhybrid.results=RNAhybrid.results[order(RNAhybrid.results$energy,decreasing = F),]



#which(RIsearch.results$TP == TRUE)

RIsearch.ranks=unlist(mapply(sRNA_TP_rank_find,sRNA_locus_tags,MoreArgs = list(tool_DF=RIsearch.results),USE.NAMES = F))
IntaRNA.ranks=unlist(mapply(sRNA_TP_rank_find,sRNA_locus_tags,MoreArgs = list(tool_DF=IntaRNA.results),USE.NAMES = F))
RNAplex.ranks=unlist(mapply(sRNA_TP_rank_find,sRNA_locus_tags,MoreArgs = list(tool_DF=RNAplex.results),USE.NAMES = F))
RNAcofold.ranks=unlist(mapply(sRNA_TP_rank_find,sRNA_locus_tags,MoreArgs = list(tool_DF=RNAcofold.results),USE.NAMES = F))
Pairfold.ranks=unlist(mapply(sRNA_TP_rank_find,sRNA_locus_tags,MoreArgs = list(tool_DF=Pairfold.results),USE.NAMES = F))
RNAup.ranks=unlist(mapply(sRNA_TP_rank_find,sRNA_locus_tags,MoreArgs = list(tool_DF=RNAup.results),USE.NAMES = F))
RNAduplex.ranks=unlist(mapply(sRNA_TP_rank_find,sRNA_locus_tags,MoreArgs = list(tool_DF=RNAduplex.results),USE.NAMES = F))
RNAhybrid.ranks=unlist(mapply(sRNA_TP_rank_find,sRNA_locus_tags,MoreArgs = list(tool_DF=RNAhybrid.results),USE.NAMES = F))




#boxplot(IntaRNA.ranks,RNAup.ranks,RIsearch.ranks,RNAplex.ranks,RNAcofold.ranks,Pairfold.ranks,names=c("IntaRNA","RNAup","RIsearch","RNAplex","RNAcofold","Pairfold"),ylab="Ranking")
#title("RNA-RNA interaction tools")


RIsearch.plot.numbers=mapply(target_per_sRNA,0:100,list(tool.result.df=RIsearch.results))
IntaRNA.plot.numbers=mapply(target_per_sRNA,0:100,list(tool.result.df=IntaRNA.results))
RNAplex.plot.numbers=mapply(target_per_sRNA,0:100,list(tool.result.df=RNAplex.results))
RNAcofold.plot.numbers=mapply(target_per_sRNA,0:100,list(tool.result.df=RNAcofold.results))
Pairfold.plot.numbers=mapply(target_per_sRNA,0:100,list(tool.result.df=Pairfold.results))
RNAup.plot.numbers=mapply(target_per_sRNA,0:100,list(tool.result.df=RNAup.results))
RNAduplex.plot.numbers=mapply(target_per_sRNA,0:100,list(tool.result.df=RNAduplex.results))
RNAhybrid.plot.numbers=mapply(target_per_sRNA,0:100,list(tool.result.df=RNAhybrid.results))


plot.new()
par(lwd=3)
plot(0:100,RIsearch.plot.numbers,type="n",ylim=c(0,60),axes=F,xlab="",ylab="")
lines(RIsearch.plot.numbers,col="red")
lines(IntaRNA.plot.numbers,col="blue")
lines(RNAplex.plot.numbers,col="darkgreen")
lines(RNAcofold.plot.numbers,col="black")
lines(Pairfold.plot.numbers,col="purple")
lines(RNAup.plot.numbers,col="orange")
lines(RNAduplex.plot.numbers,col="brown")
lines(RNAhybrid.plot.numbers,col="dimgrey")
axis(1)
axis(2)
mtext("# Target predictions per sRNA",side=1,line=2.5)
mtext("# True positives",side=2,line=2.5)
legend("topleft",legend=c("IntaRNA","RNAup","RNAplex","RNAduplex","RIsearch","RNAhybrid","Pairfold","RNAcofold"),fill=c("blue","orange","darkgreen","brown","red","dimgrey","purple","black"),bty="n")











