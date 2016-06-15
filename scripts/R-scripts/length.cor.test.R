#!/usr/bin/Rscript
arguments=commandArgs(TRUE)
require(evir)
require(evd)
require("MASS")

mcc.corr=read.csv(arguments[1],header=F,sep="\t")

cat(cor.test(mcc.corr$V1,mcc.corr$V2,method="pearson")$estimate,"\t",cor.test(mcc.corr$V1,mcc.corr$V2,method="pearson",alternative = "less")$p.value,"\n")