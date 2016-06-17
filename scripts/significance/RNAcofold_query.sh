#!/bin/bash

#$ -N RNAcofold
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash

# $1 sRNA file
# $2 Target file




program="RNAcofold"


esl-shuffle -N 200 -d --seed 255 $2 > $1.$2.$program.200Shuffled.fasta

cat $2 $1.$2.$program.200Shuffled.fasta > $1.$2.$program.Native.200Shuffled.fasta


RNAxfold_parser.py -s $1 -target $1.$2.$program.Native.200Shuffled.fasta -c $1.$2.interactions.fasta;






/home/suu13/progs/ViennaRNA-2.1.9/Progs/RNAcofold --noPS < $1.$2.interactions.fasta > $1.$2.RNAcofold.result.txt;
gawk '{if(match($0,/ \((.*)\)/,m)) print m[1]; else if(/>/) printf $1"\t";}' $1.$2.RNAcofold.result.txt | awk '{sub(/>/,"",$0); sub(/-/,"\t",$0); print $0}' > $1.$2.RNAcofold.result.csv;


