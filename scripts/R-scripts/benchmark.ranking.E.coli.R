require("RColorBrewer")
require("pracma")

Colors=brewer.pal(12, "Paired")

roc_like=function(x,tool.result.df){
  
  return(as.numeric(table(tool.result.df$V3<=x)['TRUE']))
  
}

significant_count=function(df){
  
  return(c(as.numeric(table(df$V1 <0.05)['TRUE'])/length(df$V1),as.numeric(table(df$V2 <0.05)['TRUE'])/length(df$V2)))
  
}


RIsearch.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/RIsearch/RIsearch.results.csv")
IntaRNA.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/IntaRNA/IntaRNA.results.csv")
RNAplex.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/RNAplex/RNAplex.results.csv")
RNAcofold.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/RNAcofold/RNAcofold.results.csv")
Pairfold.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/pairfold/pairfold.results.csv")
RNAup.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/RNAup/RNAup.results.csv")
RNAduplex.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/RNAduplex/RNAduplex.results.csv")
RNAhybrid.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/RNAhybrid/RNAhybrid.results.csv")
bifold.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/bifold/bifold.results.csv")
DuplexFold.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/DuplexFold/DuplexFold.results.csv")
ssearch.results=read.table("/home/suu13/projects/benchmark/bacteria/selected/ssearch36/ssearch36.results.csv")

IntaRNA.results.1000=read.table("/home/suu13/projects/benchmark/bacteria/selected/IntaRNA/IntaRNA.results1000.csv")
RIsearch.results.1000=read.table("/home/suu13/projects/benchmark/bacteria/selected/RIsearch/RIsearch.results1000.csv")






significant_results_matrix=matrix(c(significant_count(RIsearch.results),
                                significant_count(IntaRNA.results),
                                significant_count(RNAplex.results),
                                significant_count(RNAcofold.results),
                                significant_count(RNAup.results),
                                significant_count(RNAduplex.results),
                                significant_count(RNAhybrid.results),
                                significant_count(Pairfold.results),
                                significant_count(DuplexFold.results),
                                significant_count(bifold.results),
                                significant_count(ssearch.results)
                                ),
                              nrow=2)
#rownames(significant_results_matrix)=c("Extreme","Normal")



RIsearch.plot.numbers=mapply(roc_like,1:100,list(tool.result.df=RIsearch.results))
IntaRNA.plot.numbers=mapply(roc_like,1:100,list(tool.result.df=IntaRNA.results))
RNAplex.plot.numbers=mapply(roc_like,1:100,list(tool.result.df=RNAplex.results))
RNAcofold.plot.numbers=mapply(roc_like,1:100,list(tool.result.df=RNAcofold.results))
Pairfold.plot.numbers=mapply(roc_like,1:100,list(tool.result.df=Pairfold.results))
RNAup.plot.numbers=mapply(roc_like,1:100,list(tool.result.df=RNAup.results))
RNAduplex.plot.numbers=mapply(roc_like,1:100,list(tool.result.df=RNAduplex.results))
RNAhybrid.plot.numbers=mapply(roc_like,1:100,list(tool.result.df=RNAhybrid.results))
bifold.plot.numbers=mapply(roc_like,1:100,list(tool.result.df=bifold.results))
DuplexFold.plot.numbers=mapply(roc_like,1:100,list(tool.result.df=DuplexFold.results))
ssearch.plot.numbers=mapply(roc_like,1:100,list(tool.result.df=ssearch.results))



RIsearch.plot.AUC=trapz(1:100,RIsearch.plot.numbers)
IntaRNA.plot.AUC=trapz(1:100,IntaRNA.plot.numbers)
RNAplex.plot.AUC=trapz(1:100,RNAplex.plot.numbers)
RNAcofold.plot.AUC=trapz(1:100,RNAcofold.plot.numbers)
Pairfold.plot.AUC=trapz(1:100,Pairfold.plot.numbers)
RNAup.plot.AUC=trapz(1:100,RNAup.plot.numbers)
RNAduplex.plot.AUC=trapz(1:100,RNAduplex.plot.numbers)
RNAhybrid.plot.AUC=trapz(1:100,RNAhybrid.plot.numbers)
bifold.plot.AUC=trapz(1:100,bifold.plot.numbers)
DuplexFold.plot.AUC=trapz(1:100,DuplexFold.plot.numbers)
ssearch.plot.AUC=trapz(1:100,ssearch.plot.numbers)

final_table=data.frame(AUC=c(RIsearch.plot.AUC,
                             IntaRNA.plot.AUC,
                             RNAplex.plot.AUC,
                             RNAcofold.plot.AUC,
                             Pairfold.plot.AUC,
                             RNAup.plot.AUC,
                             RNAduplex.plot.AUC,
                             RNAhybrid.plot.AUC,
                             bifold.plot.AUC,
                             DuplexFold.plot.AUC,
                             ssearch.plot.AUC))

final_table=signif(cbind(final_table,t(significant_results_matrix)),digits = 2)

final_table=data.frame(AUC=final_table[,1],TPR=paste0(final_table[,2],"(",final_table[,3],")"))
rownames(final_table)=c("RIsearch","IntaRNA","RNAplex","RNAcofold","Pairfold","RNAup","RNAduplex","RNAhybrid","bifold","DuplexFold","ssearch")



boxplot(data.frame(RIsearch=RIsearch.results$V3,
                   IntaRNA=IntaRNA.results$V3,
                   RNAplex=RNAplex.results$V3,
                   RNAcofold=RNAcofold.results$V3,
                   RNAup=RNAup.results$V3,
                   Pairfold=Pairfold.results$V3,
                   RNAduplex=RNAduplex.results$V3,
                   RNAhybrid=RNAhybrid.results$V3,
                   bifold=bifold.results$V3,
                   DuplexFold=DuplexFold.results$V3,
                   ssearch=ssearch.results$V3))
#barplot
barplot(significant_results_matrix,beside = T,ylim=c(0,50))

#roc like plot


pdf("bacterial.rankings.pdf")
par(lwd=4)
plot(1:100,RIsearch.plot.numbers,type="n",ylim=c(1,250),axes=F,xlab="",ylab="")
lines(RIsearch.plot.numbers,col=Colors[1])
lines(IntaRNA.plot.numbers,col=Colors[2])
lines(RNAplex.plot.numbers,col=Colors[3])
lines(RNAcofold.plot.numbers,col=Colors[4])
lines(Pairfold.plot.numbers,col=Colors[5])
lines(RNAup.plot.numbers,col=Colors[6])
lines(RNAduplex.plot.numbers,col=Colors[7])
lines(RNAhybrid.plot.numbers,col=Colors[8])
lines(bifold.plot.numbers,col=Colors[9])
lines(DuplexFold.plot.numbers,col=Colors[10])
lines(ssearch.plot.numbers,col=Colors[11])

axis(1)
axis(2)
mtext("# Target predictions per sRNA",side=1,line=2.5)
mtext("# True positives",side=2,line=2.5)


legend("topleft",legend=c("RIsearch","IntaRNA","RNAplex","RNAcofold","Pairfold","RNAup","RNAduplex","RNAhybrid","bifold","DuplexFold","ssearch"),
       col=Colors[1:11],bty="n",lty=1,lwd=4)


dev.off()



write.table(final_table,file="/home/suu13/projects/benchmark/bacteria/bacterial.results.table.csv",sep="\t")






