#!/bin/bash

#$ -N RNAduplex
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash

# $1 sRNA file
# $2 Target file



if [ ! -f $2.200Shuffled.fasta ]; then esl-shuffle -N 200 -d --seed 255 $2 > $2.200Shuffled.fasta; fi

if [ ! -f $2.Native.200Shuffled.fasta ]; then cat $2 $2.200Shuffled.fasta > $2.Native.200Shuffled.fasta; fi



RNAxfold_parser.py -s $1 -targetRNA $2.Native.200Shuffled.fasta -duplex $1.$2.RNAduplex.fasta;


RNAduplex --noLP < $1.$2.RNAduplex.fasta > $1.$2.RNAduplex.result.txt;

gawk '{match($0,/ \((.*)\)/,m); if(NR%3 ==0) printf m[1]"\n"; else printf $1"\t"}' $1.$2.RNAduplex.result.txt | awk '{gsub(/>/,"",$0); print $0}' > $1.$2.RNAduplex.result.csv;
