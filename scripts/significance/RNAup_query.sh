#!/bin/bash

#$ -N RNAup
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash

# $1 sRNA file
# $2 Target file
# $3 
# $4
# $5 CPU

program="RNAup"


esl-shuffle -N 200 -d --seed 255 $2 > $1.$2.$program.200Shuffled.fasta

cat $2 $1.$2.$program.200Shuffled.fasta > $1.$2.$program.Native.200Shuffled.fasta

RNA_prediction_wrapper.py -w $3 $4 -sRNA $1 -targetRNA $1.$2.$program.Native.200Shuffled.fasta -program RNAup -c $5 | sort -k2,2 -V > $1.$2.RNAup.result.parsed.csv;
