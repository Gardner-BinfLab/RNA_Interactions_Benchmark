#!/bin/bash

#$ -N pairfold
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash
#$ -l mf=15G

# $1 sRNA file
# $2 Target file




if [ ! -f $2.200Shuffled.fasta ]; then esl-shuffle -N 200 -d --seed 255 $2 > $2.200Shuffled.fasta; fi

if [ ! -f $2.Native.200Shuffled.fasta ]; then cat $2 $2.200Shuffled.fasta > $2.Native.200Shuffled.fasta; fi


RNAxfold_parser.py -s $1 -targetRNA $2.Native.200Shuffled.fasta > $1.$2.pairfold.result.csv;
