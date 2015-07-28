#!/bin/bash

#$ -N RNAup
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash

# $1 sRNA file
# $2 Target file



if [ ! -f $2.200Shuffled.fasta ]; then esl-shuffle -N 200 -d --seed 255 $2 > $2.200Shuffled.fasta; fi

cat $1 $2 $2.200Shuffled.fasta > $1.$2.Native.200Shuffled.fasta;

FASTA2ID_Fix.py -f $1.$2.Native.200Shuffled.fasta -v | sponge $1.$2.Native.200Shuffled.fasta; #convert vienna

RNAup -b -o --interaction_first < $1.$2.Native.200Shuffled.fasta > $1.$2.RNAup.result.txt;

RNAup_allparser.py $1.$2.RNAup.result.txt $1.$2.RNAup.result.parsed.csv;
