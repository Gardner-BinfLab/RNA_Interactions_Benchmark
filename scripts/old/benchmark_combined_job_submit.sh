#!/bin/bash

# $1 sRNA
# $2 target RNA

for i in /home/suu13/bin/*query.sh; do ln -s $i; done
ln -s /home/suu13/progs/fasta-36.3.6d/bin/ssearch36
ln -s /home/suu13/progs/fasta-36.3.6d/data/rna.mat

sRNA_length=$(cat $1 | grep -v ">" | awk '{print length($0)}')

if [ $sRNA_length  -gt 31 ]; then sRNA_length=31; fi; #31 is internal default of RNAplex

qsub IntaRNA_query.sh $1 $2
qsub RNAcofold_query.sh $1 $2
qsub RNAhybrid_query.sh $1 $2
qsub RIsearch_query.sh $1 $2
qsub RNAduplex_query.sh $1 $2
qsub RNAplex_query.sh $1 $2 $sRNA_length
qsub -pe multi_thread 10 RNAup_query.sh $1 $2 10
qsub -pe multi_thread 10 pairfold_query.sh $1 $2 10
qsub -pe multi_thread 2 DuplexFold_query.sh $1 $2
qsub -pe multi_thread 20 bifold_query.sh $1 $2 20
./ssearch36_query.sh $1 $2

