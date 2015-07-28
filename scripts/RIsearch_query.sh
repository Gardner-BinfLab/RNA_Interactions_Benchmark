#!/bin/bash

#$ -N RIsearch
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash

# $1 sRNA file
# $2 Target file


if [ ! -f $2.200Shuffled.fasta ]; then esl-shuffle -N 200 -d --seed 255 $2 > $2.200Shuffled.fasta; fi

cat $2 $2.200Shuffled.fasta > $2.Native.200Shuffled.fasta;


RIsearch -q $1 -t $2.Native.200Shuffled.fasta -p2 -d 30 > $1.$2.RIsearch.result.txt;




