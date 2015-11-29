require("RColorBrewer")
require("pracma")


tools=c("RIsearch","IntaRNA","RNAplex","RNAcofold","pairfold","RNAup","RNAduplex","RNAhybrid","bifold","DuplexFold","ssearch36")

significant_count=function(df){
  
  return(c(as.numeric(table(df$V1 <0.05)['TRUE'])/length(df$V1),as.numeric(table(df$V2 <0.05)['TRUE'])/length(df$V2)))
  
}

overall_significant_count=function(df){
  
  df=rbind(df[[1]],df[[2]],df[[3]])
  return(c(as.numeric(table(df$V1 <0.05)['TRUE'])/length(df$V1),as.numeric(table(df$V2 <0.05)['TRUE'])/length(df$V2)))
  
}


read_results_mirna=function(species,tool) {
  return(read.table(paste0("/home/suu13/projects/benchmark/eukaryotes/results/",species,".miRNA.",tool,".results.csv")))
  
}

read_results_sirna=function(species,tool) {
  return(read.table(paste0("/home/suu13/projects/benchmark/eukaryotes/results/",species,".siRNA.",tool,".results.csv")))
  
}

read_results_pirna=function(species,tool) {
  return(read.table(paste0("/home/suu13/projects/benchmark/eukaryotes/results/",species,".piRNA.",tool,".results.csv")))
  
}

species_list_mirna=c("elegans","human","arabidopsis")


RIsearch.results=lapply(species_list_mirna,read_results_mirna,tool="RIsearch")
IntaRNA.results=lapply(species_list_mirna,read_results_mirna,tool="IntaRNA")
RNAplex.results=lapply(species_list_mirna,read_results_mirna,tool="RNAplex")
RNAcofold.results=lapply(species_list_mirna,read_results_mirna,tool="RNAcofold")
Pairfold.results=lapply(species_list_mirna,read_results_mirna,tool="pairfold")
RNAup.results=lapply(species_list_mirna,read_results_mirna,tool="RNAup")
RNAduplex.results=lapply(species_list_mirna,read_results_mirna,tool="RNAduplex")
RNAhybrid.results=lapply(species_list_mirna,read_results_mirna,tool="RNAhybrid")
bifold.results=lapply(species_list_mirna,read_results_mirna,tool="bifold")
DuplexFold.results=lapply(species_list_mirna,read_results_mirna,tool="DuplexFold")
ssearch.results=lapply(species_list_mirna,read_results_mirna,tool="ssearch36")



mirna_results_df=data.frame(rbind(
c(significant_count(RIsearch.results[[1]]),significant_count(RIsearch.results[[2]]),significant_count(RIsearch.results[[3]]),significant_count(read_results_sirna("arabidopsis","RIsearch")),significant_count(read_results_pirna("mouse","RIsearch")),overall_significant_count(RIsearch.results)),
c(significant_count(IntaRNA.results[[1]]),significant_count(IntaRNA.results[[2]]),significant_count(IntaRNA.results[[3]]),significant_count(read_results_sirna("arabidopsis","IntaRNA")),significant_count(read_results_pirna("mouse","IntaRNA")),overall_significant_count(IntaRNA.results)),
c(significant_count(RNAplex.results[[1]]),significant_count(RNAplex.results[[2]]),significant_count(RNAplex.results[[3]]),significant_count(read_results_sirna("arabidopsis","RNAplex")),significant_count(read_results_pirna("mouse","RNAplex")),overall_significant_count(RNAplex.results)),
c(significant_count(RNAcofold.results[[1]]),significant_count(RNAcofold.results[[2]]),significant_count(RNAcofold.results[[3]]),significant_count(read_results_sirna("arabidopsis","RNAcofold")),significant_count(read_results_pirna("mouse","RNAcofold")),overall_significant_count(RNAcofold.results)),
c(significant_count(Pairfold.results[[1]]),significant_count(Pairfold.results[[2]]),significant_count(Pairfold.results[[3]]),significant_count(read_results_sirna("arabidopsis","pairfold")),significant_count(read_results_pirna("mouse","pairfold")),overall_significant_count(Pairfold.results)),
c(significant_count(RNAup.results[[1]]),significant_count(RNAup.results[[2]]),significant_count(RNAup.results[[3]]),significant_count(read_results_sirna("arabidopsis","RNAup")),significant_count(read_results_pirna("mouse","RNAup")),overall_significant_count(RNAup.results)),
c(significant_count(RNAduplex.results[[1]]),significant_count(RNAduplex.results[[2]]),significant_count(RNAduplex.results[[3]]),significant_count(read_results_sirna("arabidopsis","RNAduplex")),significant_count(read_results_pirna("mouse","RNAduplex")),overall_significant_count(RNAduplex.results)),
c(significant_count(RNAhybrid.results[[1]]),significant_count(RNAhybrid.results[[2]]),significant_count(RNAhybrid.results[[3]]),significant_count(read_results_sirna("arabidopsis","RNAhybrid")),significant_count(read_results_pirna("mouse","RNAhybrid")),overall_significant_count(RNAhybrid.results)),
c(significant_count(bifold.results[[1]]),significant_count(bifold.results[[2]]),significant_count(bifold.results[[3]]),significant_count(read_results_sirna("arabidopsis","bifold")),significant_count(read_results_pirna("mouse","bifold")),overall_significant_count(bifold.results)),
c(significant_count(DuplexFold.results[[1]]),significant_count(DuplexFold.results[[2]]),significant_count(DuplexFold.results[[3]]),significant_count(read_results_sirna("arabidopsis","DuplexFold")),significant_count(read_results_pirna("mouse","DuplexFold")),overall_significant_count(DuplexFold.results)),
c(significant_count(ssearch.results[[1]]),significant_count(ssearch.results[[2]]),significant_count(ssearch.results[[3]]),significant_count(read_results_sirna("arabidopsis","ssearch36")),significant_count(read_results_pirna("mouse","ssearch36")),overall_significant_count(ssearch.results))
))

mirna_results_df[is.na(mirna_results_df)]=0





formatted_mirna_results_df=signif(mirna_results_df,digits = 2)
formatted_mirna_results_df=data.frame(Elegans=paste0(formatted_mirna_results_df[,1],"(",formatted_mirna_results_df[,2],")"),
                                      Human=paste0(formatted_mirna_results_df[,3],"(",formatted_mirna_results_df[,4],")"),
                                      Arabidopsis=paste0(formatted_mirna_results_df[,5],"(",formatted_mirna_results_df[,6],")"),
                                      Arabidopsis_siRNA=paste0(formatted_mirna_results_df[,7],"(",formatted_mirna_results_df[,8],")"),
                                      Mouse_piRNA=paste0(formatted_mirna_results_df[,9],"(",formatted_mirna_results_df[,10],")"),
                                      Overall=paste0(formatted_mirna_results_df[,11],"(",formatted_mirna_results_df[,12],")"))

rownames(formatted_mirna_results_df)=c("RIsearch","IntaRNA","RNAplex","RNAcofold","Pairfold","RNAup","RNAduplex","RNAhybrid","bifold","DuplexFold","ssearch")


write.table(formatted_mirna_results_df,file="/home/suu13/projects/benchmark/eukaryotes/results/eukaryotic.mirna.results.table.csv",sep="\t")

