require("RColorBrewer")
require("pracma")


tools=c("RIsearch","IntaRNA","RNAplex","RNAcofold","pairfold","RNAup","RNAduplex","RNAhybrid","bifold","DuplexFold","ssearch36")

TPR_calculate_full=function(df){
  
  return(c(as.numeric(table(df$P.Gumb <0.05)['TRUE'])/length(df$P.Gumb),as.numeric(table(df$P.Norm <0.05)['TRUE'])/length(P.Norm)))
  
}

TPR_calculate=function(df){
  
  return(as.numeric(table(df$P.Norm <0.05)['TRUE'])/length(df$P.Norm))
  
}

PPV_calculate=function(result_list){
  
  df=do.call(rbind,result_list)
  return(c(as.numeric(table(df$P.Norm <0.05)['TRUE'])/sum(df$Norm.FP)))
  
}

overall_TPR_calculate=function(result_list){
  
  df=do.call(rbind,result_list)
  return(c(as.numeric(table(df$P.Norm <0.05)['TRUE'])/length(df$P.Norm)))
  
}


read_results_generic=function(species,tool) {
  df=read.table(paste0("/home/suu13/projects/benchmark/eukaryotes/results/",species,tool,".results.csv"))
  colnames(df)=c("P.Gumb","P.Norm","Rank","Z.score","Gum.FP","Norm.FP")
  return(df)
  
}


species_list_mirna=c("elegans.miRNA.","human.miRNA.","arabidopsis.miRNA.")
species_list_snRNA=c("elegans.snRNA.","human.snRNA.","arabidopsis.snRNA.","yeast.snRNA.")
species_list_pirna=c("mouse.piRNA.")
species_list_snoRNA=c("elegans.snoRNA.","human.snoRNA.","arabidopsis.snoRNA.","yeast.snoRNA.")

RIsearch.results=c(sapply(species_list_mirna,read_results_generic,tool="RIsearch",USE.NAMES=T,simplify = F),
                   sapply(species_list_snRNA,read_results_generic,tool="RIsearch",USE.NAMES=T,simplify = F),
                    sapply(species_list_pirna,read_results_generic,tool="RIsearch",USE.NAMES=T,simplify = F),  
                    sapply(species_list_snoRNA,read_results_generic,tool="RIsearch",USE.NAMES=T,simplify = F))

IntaRNA.results=c(sapply(species_list_mirna,read_results_generic,tool="IntaRNA",USE.NAMES=T,simplify = F),
                   sapply(species_list_snRNA,read_results_generic,tool="IntaRNA",USE.NAMES=T,simplify = F),
                   sapply(species_list_pirna,read_results_generic,tool="IntaRNA",USE.NAMES=T,simplify = F),  
                   sapply(species_list_snoRNA,read_results_generic,tool="IntaRNA",USE.NAMES=T,simplify = F))

RNAplex.results=c(sapply(species_list_mirna,read_results_generic,tool="RNAplex",USE.NAMES=T,simplify = F),
                  sapply(species_list_snRNA,read_results_generic,tool="RNAplex",USE.NAMES=T,simplify = F),
                  sapply(species_list_pirna,read_results_generic,tool="RNAplex",USE.NAMES=T,simplify = F),  
                  sapply(species_list_snoRNA,read_results_generic,tool="RNAplex",USE.NAMES=T,simplify = F))

RNAcofold.results=c(sapply(species_list_mirna,read_results_generic,tool="RNAcofold",USE.NAMES=T,simplify = F),
                  sapply(species_list_snRNA,read_results_generic,tool="RNAcofold",USE.NAMES=T,simplify = F),
                  sapply(species_list_pirna,read_results_generic,tool="RNAcofold",USE.NAMES=T,simplify = F),  
                  sapply(species_list_snoRNA,read_results_generic,tool="RNAcofold",USE.NAMES=T,simplify = F))

Pairfold.results=c(sapply(species_list_mirna,read_results_generic,tool="pairfold",USE.NAMES=T,simplify = F),
                    sapply(species_list_snRNA,read_results_generic,tool="pairfold",USE.NAMES=T,simplify = F),
                    sapply(species_list_pirna,read_results_generic,tool="pairfold",USE.NAMES=T,simplify = F),  
                    sapply(species_list_snoRNA,read_results_generic,tool="pairfold",USE.NAMES=T,simplify = F))


RNAup.results=c(sapply(species_list_mirna,read_results_generic,tool="RNAup",USE.NAMES=T,simplify = F),
                   sapply(species_list_snRNA,read_results_generic,tool="RNAup",USE.NAMES=T,simplify = F),
                   sapply(species_list_pirna,read_results_generic,tool="RNAup",USE.NAMES=T,simplify = F),  
                   sapply(species_list_snoRNA,read_results_generic,tool="RNAup",USE.NAMES=T,simplify = F))

RNAduplex.results=c(sapply(species_list_mirna,read_results_generic,tool="RNAduplex",USE.NAMES=T,simplify = F),
                sapply(species_list_snRNA,read_results_generic,tool="RNAduplex",USE.NAMES=T,simplify = F),
                sapply(species_list_pirna,read_results_generic,tool="RNAduplex",USE.NAMES=T,simplify = F),  
                sapply(species_list_snoRNA,read_results_generic,tool="RNAduplex",USE.NAMES=T,simplify = F))

RNAhybrid.results=c(sapply(species_list_mirna,read_results_generic,tool="RNAhybrid",USE.NAMES=T,simplify = F),
                    sapply(species_list_snRNA,read_results_generic,tool="RNAhybrid",USE.NAMES=T,simplify = F),
                    sapply(species_list_pirna,read_results_generic,tool="RNAhybrid",USE.NAMES=T,simplify = F),  
                    sapply(species_list_snoRNA,read_results_generic,tool="RNAhybrid",USE.NAMES=T,simplify = F))

bifold.results=c(sapply(species_list_mirna,read_results_generic,tool="bifold",USE.NAMES=T,simplify = F),
                    sapply(species_list_snRNA,read_results_generic,tool="bifold",USE.NAMES=T,simplify = F),
                    sapply(species_list_pirna,read_results_generic,tool="bifold",USE.NAMES=T,simplify = F),  
                    sapply(species_list_snoRNA,read_results_generic,tool="bifold",USE.NAMES=T,simplify = F))

DuplexFold.results=c(sapply(species_list_mirna,read_results_generic,tool="DuplexFold",USE.NAMES=T,simplify = F),
                 sapply(species_list_snRNA,read_results_generic,tool="DuplexFold",USE.NAMES=T,simplify = F),
                 sapply(species_list_pirna,read_results_generic,tool="DuplexFold",USE.NAMES=T,simplify = F),  
                 sapply(species_list_snoRNA,read_results_generic,tool="DuplexFold",USE.NAMES=T,simplify = F))

ssearch.results=c(sapply(species_list_mirna,read_results_generic,tool="ssearch36",USE.NAMES=T,simplify = F),
                     sapply(species_list_snRNA,read_results_generic,tool="ssearch36",USE.NAMES=T,simplify = F),
                     sapply(species_list_pirna,read_results_generic,tool="ssearch36",USE.NAMES=T,simplify = F),  
                     sapply(species_list_snoRNA,read_results_generic,tool="ssearch36",USE.NAMES=T,simplify = F))



#mirna_results_df=data.frame(rbind(
#c(TPR_calculate(RIsearch.results[[1]]),TPR_calculate(RIsearch.results[[2]]),TPR_calculate(RIsearch.results[[3]]),TPR_calculate(read_results_generic("arabidopsis.siRNA.","RIsearch")),TPR_calculate(read_results_generic("mouse.piRNA.","RIsearch")),overall_TPR_calculate(RIsearch.results)),
#c(TPR_calculate(IntaRNA.results[[1]]),TPR_calculate(IntaRNA.results[[2]]),TPR_calculate(IntaRNA.results[[3]]),TPR_calculate(read_results_generic("arabidopsis.siRNA.","IntaRNA")),TPR_calculate(read_results_generic("mouse.piRNA.","IntaRNA")),overall_TPR_calculate(IntaRNA.results)),
#c(TPR_calculate(RNAplex.results[[1]]),TPR_calculate(RNAplex.results[[2]]),TPR_calculate(RNAplex.results[[3]]),TPR_calculate(read_results_generic("arabidopsis.siRNA.","RNAplex")),TPR_calculate(read_results_generic("mouse.piRNA.","RNAplex")),overall_TPR_calculate(RNAplex.results)),
#c(TPR_calculate(RNAcofold.results[[1]]),TPR_calculate(RNAcofold.results[[2]]),TPR_calculate(RNAcofold.results[[3]]),TPR_calculate(read_results_generic("arabidopsis.siRNA.","RNAcofold")),TPR_calculate(read_results_generic("mouse.piRNA.","RNAcofold")),overall_TPR_calculate(RNAcofold.results)),
#c(TPR_calculate(Pairfold.results[[1]]),TPR_calculate(Pairfold.results[[2]]),TPR_calculate(Pairfold.results[[3]]),TPR_calculate(read_results_generic("arabidopsis.siRNA.","pairfold")),TPR_calculate(read_results_generic("mouse.piRNA.","pairfold")),overall_TPR_calculate(Pairfold.results)),
#c(TPR_calculate(RNAup.results[[1]]),TPR_calculate(RNAup.results[[2]]),TPR_calculate(RNAup.results[[3]]),TPR_calculate(read_results_generic("arabidopsis.siRNA.","RNAup")),TPR_calculate(read_results_generic("mouse.piRNA.","RNAup")),overall_TPR_calculate(RNAup.results)),
#(TPR_calculate(RNAduplex.results[[1]]),TPR_calculate(RNAduplex.results[[2]]),TPR_calculate(RNAduplex.results[[3]]),TPR_calculate(read_results_generic("arabidopsis.siRNA.","RNAduplex")),TPR_calculate(read_results_generic("mouse.piRNA.","RNAduplex")),overall_TPR_calculate(RNAduplex.results)),
#c(TPR_calculate(RNAhybrid.results[[1]]),TPR_calculate(RNAhybrid.results[[2]]),TPR_calculate(RNAhybrid.results[[3]]),TPR_calculate(read_results_generic("arabidopsis.siRNA.","RNAhybrid")),TPR_calculate(read_results_generic("mouse.piRNA.","RNAhybrid")),overall_TPR_calculate(RNAhybrid.results)),
#c(TPR_calculate(bifold.results[[1]]),TPR_calculate(bifold.results[[2]]),TPR_calculate(bifold.results[[3]]),TPR_calculate(read_results_generic("arabidopsis.siRNA.","bifold")),TPR_calculate(read_results_generic("mouse.piRNA.","bifold")),overall_TPR_calculate(bifold.results)),
#c(TPR_calculate(DuplexFold.results[[1]]),TPR_calculate(DuplexFold.results[[2]]),TPR_calculate(DuplexFold.results[[3]]),TPR_calculate(read_results_generic("arabidopsis.siRNA.","DuplexFold")),TPR_calculate(read_results_generic("mouse.piRNA.","DuplexFold")),overall_TPR_calculate(DuplexFold.results)),
#c(TPR_calculate(ssearch.results[[1]]),TPR_calculate(ssearch.results[[2]]),TPR_calculate(ssearch.results[[3]]),TPR_calculate(read_results_generic("arabidopsis.siRNA.","ssearch36")),TPR_calculate(read_results_generic("mouse.piRNA.","ssearch36")),overall_TPR_calculate(ssearch.results))))

#1 mirna elegans
#2 mirna human
#3 mirna arabidopsis
#4 snrna elegans
#5 snrna human
#6 snrna arabidopsis
#7 snrna yeast



mirna_results_df=data.frame(rbind(
c(unlist(lapply(RIsearch.results,TPR_calculate)),overall_TPR_calculate(RIsearch.results)),
c(unlist(lapply(IntaRNA.results,TPR_calculate)),overall_TPR_calculate(IntaRNA.results)),
c(unlist(lapply(RNAplex.results,TPR_calculate)),overall_TPR_calculate(RNAplex.results)),
c(unlist(lapply(RNAcofold.results,TPR_calculate)),overall_TPR_calculate(RNAcofold.results)),
c(unlist(lapply(Pairfold.results,TPR_calculate)),overall_TPR_calculate(Pairfold.results)),
c(unlist(lapply(RNAup.results,TPR_calculate)),overall_TPR_calculate(RNAup.results)),
c(unlist(lapply(RNAduplex.results,TPR_calculate)),overall_TPR_calculate(RNAduplex.results)),
c(unlist(lapply(RNAhybrid.results,TPR_calculate)),overall_TPR_calculate(RNAhybrid.results)),
c(unlist(lapply(bifold.results,TPR_calculate)),overall_TPR_calculate(bifold.results)),
c(unlist(lapply(DuplexFold.results,TPR_calculate)),overall_TPR_calculate(DuplexFold.results)),
c(unlist(lapply(ssearch.results,TPR_calculate)),overall_TPR_calculate(ssearch.results))))

rownames(mirna_results_df)=c("RIsearch","IntaRNA","RNAplex","RNAcofold","Pairfold","RNAup","RNAduplex","RNAhybrid","bifold","DuplexFold","ssearch")


mirna_results_df[is.na(mirna_results_df)]=0





#formatted_mirna_results_df=data.frame(Elegans=paste0(formatted_mirna_results_df[,1],"(",formatted_mirna_results_df[,2],")"),
#                                      Human=paste0(formatted_mirna_results_df[,3],"(",formatted_mirna_results_df[,4],")"),
#                                      Arabidopsis=paste0(formatted_mirna_results_df[,5],"(",formatted_mirna_results_df[,6],")"),
#                                      Arabidopsis_siRNA=paste0(formatted_mirna_results_df[,7],"(",formatted_mirna_results_df[,8],")"),
#                                      Mouse=paste0(formatted_mirna_results_df[,9],"(",formatted_mirna_results_df[,10],")"),
#                                      Overall=paste0(formatted_mirna_results_df[,11],"(",formatted_mirna_results_df[,12],")"))

write.table(signif(mirna_results_df),file="/home/suu13/projects/benchmark/eukaryotes/results/eukaryotic.interactions.results.table.csv",sep="\t")

