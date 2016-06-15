#!/bin/bash

#$ -N Benchmark
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash

# $1 sRNA file
# $2 Target file
# $3 window start
# $4 window stop


for i in RIsearch IntaRNA RNAplex RNAcofold RNAup RNAduplex RNAhybrid pairfold DuplexFold bifold ssearch ractip bistarna AccessFold NUPACK;

#for i in NUPACK;


#do 
#RNA_interaction_wrapper.py -s $1 -t $2 -w $3 $4 -program $i;

#done | grep -v Algorithm >> $1.$2.output.csv;

#for i in ractip bistarna;

#for i in RIsearch IntaRNA RNAplex RNAcofold RNAup RNAduplex RNAhybrid pairfold DuplexFold bifold ssearch ractip bistarna AccessFold;

do 
{ /usr/bin/time --quiet -f "$i\t%e" RNA_interaction_wrapper.py -s $1 -t $2 -w $3 $4 -program $i > /dev/null;};

done;
