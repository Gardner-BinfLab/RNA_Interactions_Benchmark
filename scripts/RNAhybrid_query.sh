#!/bin/bash

#$ -N RNAhyrid
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash

# $1 sRNA file
# $2 Target file



if [ ! -f $2.200Shuffled.fasta ]; then esl-shuffle -N 200 -d --seed 255 $2 > $2.200Shuffled.fasta; fi

cat $2 $2.200Shuffled.fasta > $2.Native.200Shuffled.fasta;


RNAhybrid -t $2.Native.200Shuffled.fasta -q $1 -d 0 1 -n 1000 -m 10000 > $1.$2.RNAhybrid.result.txt;

grep -e "target:\|miRNA :\|mfe:" $1.$2.RNAhybrid.result.txt | awk '{if(/target/) t=$2; else if (/miRNA/) m=$3; else if(/mfe/) print t,m,$2;}' >  $1.$2.RNAhybrid.result.csv;
