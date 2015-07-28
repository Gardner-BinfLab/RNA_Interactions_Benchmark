#!/bin/bash

#$ -N IntaRNA
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash

# $1 sRNA file
# $2 Target file



if [ ! -f $2.200Shuffled.fasta ]; then esl-shuffle -N 200 -d --seed 255 $2 > $2.200Shuffled.fasta; fi

cat $2 $2.200Shuffled.fasta > $2.Native.200Shuffled.fasta;

IntaRNA -m $1 -t $2.Native.200Shuffled.fasta -w 32 -o > $1.$2.IntaRNA.result.txt;

grep -E ">|^energy" $1.$2.IntaRNA.result.txt | gawk '{if(NF%3 ==0) print $2; else printf $1"\t"}' | awk '{gsub(/>/,"",$0); print $0}' > $1.$2.IntaRNA.result.csv;
