#!/usr/bin/Rscript
arguments=commandArgs(TRUE)


shuffled_results=read.csv(arguments[1],header=F,sep="\t")

median_FP_gumbel=median(as.numeric(shuffled_results[,5]))
median_FP_normal=median(as.numeric(shuffled_results[,6]))

median_rank=median(as.numeric(shuffled_results[,3]))



cat(median_FP_gumbel,"\t",median_FP_normal,"\t",median_rank,"\n")



