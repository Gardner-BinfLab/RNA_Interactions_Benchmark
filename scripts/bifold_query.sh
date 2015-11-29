#!/bin/bash

#$ -N bifold
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash
#$ -l mf=15G

# $1 sRNA file
# $2 Target file
# $3 CPU

export DATAPATH=$HOME/progs/RNAstructure/data_tables


program="bifold"

esl-shuffle -N 200 -d --seed 255 $2 > $1.$2.$program.200Shuffled.fasta

cat $2 $1.$2.$program.200Shuffled.fasta > $1.$2.$program.Native.200Shuffled.fasta


RNA_prediction_wrapper.py -sRNA $1 -targetRNA $1.$2.$program.Native.200Shuffled.fasta -program bifold -c $3 | sort -k2,2 -V > $1.$2.bifold.result.csv;
