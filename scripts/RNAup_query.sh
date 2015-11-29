#!/bin/bash

#$ -N RNAup
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash

# $1 sRNA file
# $2 Target file

program="RNAup"


esl-shuffle -N 200 -d --seed 255 $2 > $1.$2.$program.200Shuffled.fasta

cat $2 $1.$2.$program.200Shuffled.fasta > $1.$2.$program.Native.200Shuffled.fasta


cat $1 $1.$2.$program.Native.200Shuffled.fasta | sponge $1.$2.$program.Native.200Shuffled.fasta; #add sRNA at the beginning or miRNA etc

FASTA2ID_Fix.py -f $1.$2.$program.Native.200Shuffled.fasta -v | sponge $1.$2.$program.Native.200Shuffled.fasta; #convert vienna

/home/suu13/progs/ViennaRNA-2.1.9/Progs/RNAup -b -o --interaction_first < $1.$2.$program.Native.200Shuffled.fasta > $1.$2.RNAup.result.txt;

RNAup_allparser.py $1.$2.RNAup.result.txt $1.$2.RNAup.result.parsed.csv;
