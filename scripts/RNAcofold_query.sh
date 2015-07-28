#!/bin/bash

#$ -N RNAcofold
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash

# $1 sRNA file
# $2 Target file


if [ ! -f $2.200Shuffled.fasta ]; then esl-shuffle -N 200 -d --seed 255 $2 > $2.200Shuffled.fasta; fi

cat $2 $2.200Shuffled.fasta > $2.Native.200Shuffled.fasta;

RNAxfold_parser.py -s $1 -target  $2.Native.200Shuffled.fasta -c $1.$2.interactions.fasta;






RNAcofold --noPS < $1.$2.interactions.fasta > $1.$2.RNAcofold.result.txt;
gawk '{if(match($0,/ \((.*)\)/,m)) print m[1]; else if(/>/) printf $1"\t"$2"\t";}' $1.$2.RNAcofold.result.txt | awk '{gsub(/>/,"",$0); print $0}' > $1.$2.RNAcofold.result.csv;


