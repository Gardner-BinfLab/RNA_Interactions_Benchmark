#!/bin/bash


for i in RIsearch IntaRNA RNAplex RNAup RNAduplex RNAhybrid DuplexFold ssearch AccessFold;


do 
for a in /home/suu13/projects/benchmark/bacteria/shuffled.results/*$i*csv; 
do p-value-fit.R $a; 

done > /home/suu13/projects/benchmark/bacteria/shuffled.$i.csv;

done

