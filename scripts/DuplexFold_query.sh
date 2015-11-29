#!/bin/bash

#$ -N DuplexFold
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash
#$ -l mf=15G

# $1 sRNA file
# $2 Target file

export DATAPATH=$HOME/progs/RNAstructure/data_tables


program="DuplexFold"


esl-shuffle -N 200 -d --seed 255 $2 > $1.$2.$program.200Shuffled.fasta

cat $2 $1.$2.$program.200Shuffled.fasta > $1.$2.$program.Native.200Shuffled.fasta


RNA_prediction_wrapper.py -sRNA $1 -targetRNA $1.$2.$program.Native.200Shuffled.fasta -program DuplexFold > $1.$2.DuplexFold.result.csv;
