#!/usr/bin/Rscript
arguments=commandArgs(TRUE)
require(evir)
require(evd)
require("MASS")


parsed_result_csv=read.csv(arguments[1],header=F,sep="\t")

#assign median to NAs
if(anyNA(parsed_result_csv$V3)) {
  median_results=median(parsed_result_csv[!is.na(parsed_result_csv$V3),]$V3)
  parsed_result_csv[is.na(parsed_result_csv$V3),]$V3=median_results
}

gumbel_parameters=gumbel(-1*as.vector(parsed_result_csv$V3[2:length(parsed_result_csv$V3)]))$par.ests
normal_parameters=fitdistr(as.vector(parsed_result_csv$V3[2:length(parsed_result_csv$V3)]),densfun = "normal")

p_value_native_gumbel=pgumbel(-1*parsed_result_csv$V3[1],unname(gumbel_parameters[2]),unname(gumbel_parameters[1]),lower.tail = F)
p_value_native_normal=pnorm(parsed_result_csv$V3[1],mean=normal_parameters$estimate[1],sd=normal_parameters$estimate[2],lower.tail = T)
  
cat(p_value_native_gumbel,"\t",p_value_native_normal,"\t",rank(parsed_result_csv$V3)[1],"\n")


