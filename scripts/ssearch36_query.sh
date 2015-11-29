#!/bin/bash

#$ -N ssearch36
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash

# $1 sRNA file
# $2 Target file


program="ssearch36"

esl-shuffle -N 200 -d --seed 255 $2 > $1.$2.$program.200Shuffled.fasta

cat $2 $1.$2.$program.200Shuffled.fasta > $1.$2.$program.Native.200Shuffled.fasta


biof_converter.py -f $1.$2.$program.Native.200Shuffled.fasta -i fasta -o fasta -reversecomplement > $2.rev_comp.fasta;

./ssearch36 -s rna.mat $1 $2.rev_comp.fasta -3 -f 16 -g 10 -E 1000 -m 8C | grep -v ^# | awk '{print $1"\t"$2"\t"(-1*$12/$4)}' | sort -k3,3n | sort -u -k2,2V > $1.$2.ssearch36.result.csv;

