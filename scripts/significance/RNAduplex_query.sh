#!/bin/bash

#$ -N RNAduplex
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash

# $1 sRNA file
# $2 Target file

program="RNAduplex"


esl-shuffle -N 200 -d --seed 255 $2 > $1.$2.$program.200Shuffled.fasta

cat $2 $1.$2.$program.200Shuffled.fasta > $1.$2.$program.Native.200Shuffled.fasta



RNAxfold_parser.py -s $1 -targetRNA $1.$2.$program.Native.200Shuffled.fasta -duplex $1.$2.RNAduplex.fasta;


/home/suu13/progs/ViennaRNA-2.1.9/Progs/RNAduplex --noLP < $1.$2.RNAduplex.fasta > $1.$2.RNAduplex.result.txt;

gawk '{match($0,/ \((.*)\)/,m); if(NR%3 ==0) printf m[1]"\n"; else printf $1"\t"}' $1.$2.RNAduplex.result.txt | awk '{gsub(/>/,"",$0); print $0}' > $1.$2.RNAduplex.result.csv;
