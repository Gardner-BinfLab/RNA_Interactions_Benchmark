#!/bin/bash

#$ -N RNAplex
#$ -o /home/suu13/job_logs/
#$ -e /home/suu13/job_logs/
#$ -cwd
#$ -S /bin/bash
#$ -l mf=50G

# $1 sRNA file
# $2 Target file



if [ ! -f $2.200Shuffled.fasta ]; then esl-shuffle -N 200 -d --seed 255 $2 > $2.200Shuffled.fasta; fi

if [ ! -f $2.Native.200Shuffled.fasta ]; then cat $2 $2.200Shuffled.fasta > $2.Native.200Shuffled.fasta; fi


if [ -n "$JOB_ID" ]; 

then
FASTA2ID_Fix.py -f $1 -v | sponge $JOB_ID.$1; #convert vienna just in case
FASTA2ID_Fix.py -f $2.Native.200Shuffled.fasta -v | sponge $JOB_ID.$2.Native.200Shuffled.fasta; #convert vienna

/home/suu13/progs/ViennaRNA-2.1.9/Progs/RNAplfold -O --plex_output < $JOB_ID.$1;
/home/suu13/progs/ViennaRNA-2.1.9/Progs/RNAplfold -O --plex_output < $JOB_ID.$2.Native.200Shuffled.fasta;
/home/suu13/progs/ViennaRNA-2.1.9/Progs/RNAplex -l 30 -q $JOB_ID.$1 -t $JOB_ID.$2.Native.200Shuffled.fasta -a ./ > $1.$2.RNAplex.result.txt;


else
FASTA2ID_Fix.py -f $1 -v | sponge $1; #convert vienna just in case
FASTA2ID_Fix.py -f $2.Native.200Shuffled.fasta -v | sponge $2.Native.200Shuffled.fasta; #convert vienna
/home/suu13/progs/ViennaRNA-2.1.9/Progs/RNAplfold -O --plex_output < $1;
/home/suu13/progs/ViennaRNA-2.1.9/Progs/RNAplfold -O --plex_output < $2.Native.200Shuffled.fasta;
/home/suu13/progs/ViennaRNA-2.1.9/Progs/RNAplex -l 30 -q $1 -t $2.Native.200Shuffled.fasta -a ./ > $1.$2.RNAplex.result.txt;

fi;




gawk '{match($0,/ \((.*) =/,m); if(NR%3 ==0) printf m[1]"\n"; else printf $1"\t"}' $1.$2.RNAplex.result.txt | awk '{gsub(/>/,"",$0); print $0}' > $1.$2.RNAplex.result.csv;




