#!/bin/bash

#$ -N IntaRNA
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash

# $1 sRNA file
# $2 Target file

program="IntaRNA"


esl-shuffle -N 200 -d --seed 255 $2 > $1.$2.$program.200Shuffled.fasta

cat $2 $1.$2.$program.200Shuffled.fasta > $1.$2.$program.Native.200Shuffled.fasta



IntaRNA -m $1 -t $1.$2.$program.Native.200Shuffled.fasta -w 32 -o > $1.$2.IntaRNA.result.txt;

grep -E ">|^energy" $1.$2.IntaRNA.result.txt | gawk '{if(NF%3 ==0) print $2; else printf $1"\t"}' | awk '{gsub(/>/,"",$0); print $0}' > $1.$2.IntaRNA.result.csv;
