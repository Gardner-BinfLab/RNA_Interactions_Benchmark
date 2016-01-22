#!/usr/bin/env python2.7


'''
Created on 04/12/2015

@author: suu13
'''


import argparse
from Bio import SeqIO
import subprocess
import tempfile
import os
from math import sqrt
from re import finditer


__author__ = 'suu13'



'''RIsearch output example
HBII-239	20	47	gb-U13369.1	2	32	2310	-35.21
'''
def RIsearch_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()

    sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='RIsearch_sRNA_',suffix='.fasta',mode="w+",delete=False)
    sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.seq)))
    sRNA_temp_file.close()

    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='RIsearch_RNA_',suffix='.fasta',mode="w+",delete=False)
    RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq)))
    RNA_temp_file.close()


    shell_command="""RIsearch -q %s -t %s -d 30 -p2 | cut -f5,6 """ %(sRNA_temp_file.name,RNA_temp_file.name)
    predicted_region=subprocess.check_output(shell_command,shell=True)
    true_base_pairs=set(range(real_region[0],real_region[1]+1))
    predicted_base_pairs=set(range(int(predicted_region.split()[0].strip()),int(predicted_region.split()[1].strip())+1))
    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))

    print_results(args.program,True_Positives,False_Positives,False_Negatives)

    os.remove(RNA_temp_file.name)
    os.remove(sRNA_temp_file.name)
    return




'''IntaRNA output example
-------------------------
OUTPUT
-------------------------
>gb-U13369.1
CGACUCUUAGCGGUGGAUCACUCGGCUCGUGCGUCGAUGAAGAACGCAGCUAGCUGCGAGAAUUAAUGUGAAUUGCAGGACACAUUGAUCAUCGACACUUCGAACGCACUUGCGGCCCCGGGUUCCUCCCGGGGCUACGCCUGUCUGAGCGUCGCUU
>HBII-239
UGUGUGUUGGAGGAUGAAAGUACGGAGUGAUCCAUCGGCUAAGUGUCUUGUCACAAUGCUGACACUCAAACUGCUGACAGCACACG

                                      5'-C   U                   GGCU     GUCGAUGAAGAACGCAGCUAGCUGCGAGAAUUAAUGUGAAUUGCAGGACACAUUGAUCAUCGACACUUCGAACGCACUUGCGGCCCCGGGUUCCUCCCGGGGCUACGCCUGUCUGAGCGUCGCUU-3'
                                          GAC CUUAGC GGUGGAUCACUC    CGUGC
                                          CUG GAAUCG CUACCUAGUGAG    GCAUG
3'-GCACACGACAGUCGUCAAACUCACAGUCGUAACACUGUU   U      G                     AAAGUAGGAGGUUGUGUGU-5'

positions(target)     : 2 -- 32
positions seed(target): 18 -- 23
positions with dangle(target): 2 -- 33
positions(ncRNA)      : 20 -- 47
positions seed(ncRNA) : 25 -- 30
positions with dangle(ncRNA): 19 -- 48
ED target need: 6.13863 kcal/mol
ED ncRNA  need: 14.3114 kcal/mol
hybrid energy : -39.8 kcal/mol

energy: -19.35 kcal/mol
'''
def IntaRNA_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()

    sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='IntaRNA_sRNA_',suffix='.fasta',mode="w+",delete=False)
    sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.seq)))
    sRNA_temp_file.close()

    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='IntaRNA_RNA_',suffix='.fasta',mode="w+",delete=False)
    RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq)))
    RNA_temp_file.close()


    shell_command="""IntaRNA -m %s -t %s -w 32 -o | grep "positions(target)" | gawk 'match($0,/.*:\s+([0-9]+)\s+--\s+([0-9]+).*/,m){print m[1],m[2]}' """ %(sRNA_temp_file.name,RNA_temp_file.name)

    predicted_region=subprocess.check_output(shell_command,shell=True)
    true_base_pairs=set(range(real_region[0],real_region[1]+1))
    predicted_base_pairs=set(range(int(predicted_region.split()[0].strip()),int(predicted_region.split()[1].strip())+1))
    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))

    print_results(args.program,True_Positives,False_Positives,False_Negatives)

    os.remove(RNA_temp_file.name)
    os.remove(sRNA_temp_file.name)
    return


'''
>gb-U13369.1
>HBII-239
.((.(((...((((((...(((.((((.&.)))))))))))))........))).)). 114,141 :   1,29  (-17.40)
'''
def RNAplex_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()

    sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='RNAplex_sRNA_',suffix='.fasta',mode="w+",delete=False)
    sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.seq)))
    sRNA_temp_file.close()

    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='RNAplex_RNA_',suffix='.fasta',mode="w+",delete=False)
    RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq)))
    RNA_temp_file.close()

    #size=abs(real_region[1]-real_region[0])

    window_size=len(sRNA.seq) if len(RNA.seq) >= len(sRNA.seq) else len(RNA.seq)
    subprocess.check_output("RNAplfold -W %d -u %d -O --plex_output < %s" % (window_size,window_size,sRNA_temp_file.name),shell=True)
    subprocess.check_output("RNAplfold -W %d -u %d -O --plex_output < %s" % (window_size,window_size,RNA_temp_file.name),shell=True)
    shell_command="""RNAplex -l %d -q %s -t %s -a ./ | gawk 'match($0,/\s+([0-9]+),([0-9]+)\s+/,m){print m[1],m[2]}' """ %(window_size,sRNA_temp_file.name,RNA_temp_file.name)


    #subprocess.check_output("RNAplfold -O --plex_output < %s" % (sRNA_temp_file.name),shell=True)
    #subprocess.check_output("RNAplfold -O --plex_output < %s" % (RNA_temp_file.name),shell=True)
    #shell_command="""RNAplex -l 31 -q %s -t %s -a ./ | gawk 'match($0,/\s+([0-9]+),([0-9]+)\s+/,m){print m[1],m[2]}' """ %(sRNA_temp_file.name,RNA_temp_file.name)


    predicted_region=subprocess.check_output(shell_command,shell=True)
    true_base_pairs=set(range(real_region[0],real_region[1]+1))
    predicted_base_pairs=set(range(int(predicted_region.split()[0].strip()),int(predicted_region.split()[1].strip())+1))

    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))

    print_results(args.program,True_Positives,False_Positives,False_Negatives)

    os.remove(RNA_temp_file.name)
    os.remove(sRNA_temp_file.name)
    return

'''
>deneme
UGUGUGUUGGAGGAUGAAAGUACGGAGUGAUCCAUCGGCUAAGUGUCUUGUCACAAUGCUGACACUCAAACUGCUGACAGCACACG&CGACUCUUAGCGGUGGAUCACUCGGCUCGUGCGUCGAUGAAGAACGCAGCUAGCUGCGAGAAUUAAUGUGAAUUGCAGGACACAUUGAUCAUCGACACUUCGAACGCACUUGCGGCCCCGGGUUCCUCCCGGGGCUACGCCUGUCUGAGCGUCGCUU
.(((((((.((((......((((((((((((((((((.(((((.(((.(((((......))))).......((((...))))....&.))).))))))))))))))))))....)))))((((((((....(((((....)))))....((((((((..........)))))))))))))))).)))).)))))))..(((((((((((.....)))))))).(((((.....).))))))).. (-105.70)
'''
def RNAcofold_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()

    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='RNAcofold_RNA_',suffix='.fasta',mode="w+",delete=False)
    RNA_temp_file.write(">%s-%s\n%s&%s\n" %(str(sRNA.description),str(RNA.description),str(sRNA.seq),str(RNA.seq)))
    RNA_temp_file.close()

    shell_command="""RNAcofold --noPS < %s | grep ")" | cut -d'&' -f2 | awk '{print $1}' """ %(RNA_temp_file.name)

    predicted_region=subprocess.check_output(shell_command,shell=True)



    true_base_pairs=set(range(real_region[0],real_region[1]+1))
    predicted_base_pairs=set(Parse_Dot_Bracket(predicted_region))

    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))

    if True_Positives != 0 or False_Positives != 0:
        print_results(args.program,True_Positives,False_Positives,False_Negatives)
    else:
        print_results(args.program,True_Positives,1,False_Negatives) #add pseudo count for RNAcofold

    os.remove(RNA_temp_file.name)
    return


'''
Seq: TGTGTGTTGGAGGATGAAAGTACGGAGTGATCCATCGGCTAAGTGTCTTGTCACAATGCTGACACTCAAACTGCTGACAGCACACG CGACTCTTAGCGGTGGATCACTCGGCTCGTGCGTCGATGAAGAACGCAGCTAGCTGCGAGAATTAATGTGAATTGCAGGACACATTGATCATCGACACTTCGAACGCACTTGCGGCCCCGGGTTCCTCCCGGGGCTACGCCTGTCTGAGCGTCGCTT
MFE: .(((((((.((((......((((((((((((((((((.(((((.(((.........(((((.(((.......).)).))))).... .))).))))))))))))))))))....)))))((((((((....(((((....)))))....((((((((..........)))))))))))))))).)))).)))))))..(((((((((((.....)))))))).(((((.....).)))))))..  -110.76
'''
def Pairfold_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()


    shell_command="""/home/suu13/misc_stuff/MultiRNAFold-2.0/pairfold "%s" "%s" -m RNA | grep MFE | awk '{print $3}' """ % (str(sRNA.seq),str(RNA.seq)) #enter pairfold actual path

    predicted_region=subprocess.check_output(shell_command,shell=True)

    true_base_pairs=set(range(real_region[0],real_region[1]+1))
    predicted_base_pairs=set(Parse_Dot_Bracket(predicted_region))

    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))

    if True_Positives != 0 or False_Positives != 0:
        print_results(args.program,True_Positives,False_Positives,False_Negatives)
    else:
        print_results(args.program,True_Positives,1,False_Negatives) #add pseudo count for pairfold

    return


'''
>deneme
>deneme1
(((.((((((((((((((((((&))))))))))))).))))).)))   2,23  :  25,47  (-12.71 = -35.47 + 10.37 + 12.40)
GACUCUUAGCGGUGGAUCACUC&GAGUGAUCCAUCGGCUAAGUGUC
'''
def RNAup_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()

    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='RNAup_RNA_',suffix='.fasta',mode="w+",delete=False)
    RNA_temp_file.write(">%s\n%s\n>%s\n%s\n" %(str(sRNA.description),str(sRNA.seq),str(RNA.description),str(RNA.seq)))
    RNA_temp_file.close()

    size=abs(real_region[1]-real_region[0]) # for bacteria +5 is enough
    shell_command="""RNAup -w %d -b -o --interaction_first -3 -5 < %s | gawk 'match($0,/\s+([0-9]+),([0-9]+)\s+:\s+/,m){print m[1],m[2]}'""" %(size,RNA_temp_file.name) #-3 -5 for miRNAs and snoRNAs
    predicted_region=subprocess.check_output(shell_command,shell=True)



    true_base_pairs=set(range(real_region[0],real_region[1]+1))
    predicted_base_pairs=set(range(int(predicted_region.split()[0].strip()),int(predicted_region.split()[1].strip())+1))

    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))

    print_results(args.program,True_Positives,False_Positives,False_Negatives)

    os.remove(RNA_temp_file.name)
    return

'''
>deneme
>deneme1
.(((((((((((((.....((.(((((((.((.((.((.((((((((((((((((...((..((((....((((((((.((((.&.)))))))).........))))..)).))..))......))))....))))))))).))).)))).))))))))).))..................))))))...))).)))).   1,84  :  28,141 (-50.80)
'''
def RNAduplex_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()

    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='RNAduplex_RNA_',suffix='.fasta',mode="w+",delete=False)
    RNA_temp_file.write(">%s\n%s\n>%s\n%s\n" %(str(sRNA.description),str(sRNA.seq),str(RNA.description),str(RNA.seq)))
    RNA_temp_file.close()

    shell_command="""RNAduplex --noLP < %s | gawk 'match($0,/\s+:\s+([0-9]+),([0-9]+) /,m){print m[1],m[2]}'""" %(RNA_temp_file.name)

    predicted_region=subprocess.check_output(shell_command,shell=True)



    true_base_pairs=set(range(real_region[0],real_region[1]+1))
    predicted_base_pairs=set(range(int(predicted_region.split()[0].strip()),int(predicted_region.split()[1].strip())+1))

    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))


    print_results(args.program,True_Positives,False_Positives,False_Negatives)

    os.remove(RNA_temp_file.name)
    return


'''
target: gb-U13369.1
length: 157
miRNA : HBII-239
length: 86

mfe: -60.2 kcal/mol
p-value: 1.000000

position  28
target 5'   C         AUGAAGAAC    CUAGCUGC   AA   A       AAUU         A        CAUC          A  G    CUUGCGGCCCCGGGU      CGG   U    C 3'
             GUGC GUCG         GCAG        GAG  UUA    UGUG    GCAGGACAC    UUGAU    GA CACUUCG AC   CA               UCCUCC   GGC ACGC
             CACG CAGU         CGUC        CUC  AGU    ACAC    UGUUCUGUG    GGCUA    CU GUGAGGC UG   GU               AGGAGG   UUG UGUG
miRNA  3' GCA    A                 AAA        AC   CGUA                 AAUC     C     A       A  AAA                                  U 5'
'''
def RNAhybrid_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()

    sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='RNAhybrid_sRNA_',suffix='.fasta',mode="w+",delete=False)
    sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.seq)))
    sRNA_temp_file.close()

    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='RNAhybrid_RNA_',suffix='.fasta',mode="w+",delete=False)
    RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq)))
    RNA_temp_file.close()


    shell_command="""RNAhybrid -q %s -t %s -d 0 1 -n 10000 -m 10000 | grep "position" | awk '{print $2}' """ %(sRNA_temp_file.name,RNA_temp_file.name)


    predicted_region=subprocess.check_output(shell_command,shell=True) #returns only start of target RNA

    true_base_pairs=set(range(real_region[0],real_region[1]+1))
    predicted_base_pairs=set(range(int(predicted_region.strip()),int(predicted_region.strip())+len(sRNA.seq)+1))


    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))

    print_results(args.program,True_Positives,False_Positives,False_Negatives)

    os.remove(RNA_temp_file.name)
    os.remove(sRNA_temp_file.name)
    return


'''
  246  ENERGY = -104.4  HBII-239_gb-U13369.1
    1 T       0    2    0    1
    2 G       1    3  198    2
    3 T       2    4  197    3
    4 G       3    5  196    4
    5 T       4    6  195    5
    6 G       5    7  194    6
    7 T       6    8  193    7
    8 T       7    9  192    8
    9 G       8   10    0    9
   10 G       9   11  190   10
   11 A      10   12  189   11
   12 G      11   13  188   12
   13 G      12   14  187   13
   14 A      13   15    0   14
   15 T      14   16    0   15
   16 G      15   17    0   16
'''
def bifold_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()


    sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='bifold_sRNA_',suffix='.fasta',mode="w+",delete=False)
    sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.seq)))
    sRNA_temp_file.close()

    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='bifold_RNA_',suffix='.fasta',mode="w+",delete=False)
    RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq)))
    RNA_temp_file.close()

    shell_command="""bifold -m 1 %s %s %s > /dev/null && cat %s | grep -v ENERGY | gawk '{print $(NF-1),$NF}'""" %(sRNA_temp_file.name,RNA_temp_file.name,RNA_temp_file.name+".out",RNA_temp_file.name+".out")


    predicted_region=subprocess.check_output(shell_command,shell=True).splitlines()

    true_base_pairs=set(range(real_region[0],real_region[1]+1))

    predicted_base_pairs=list()
    for i in range(0,len(predicted_region)):
        if i+1 > len(sRNA.seq) and int(predicted_region[i].split()[0]) < len(sRNA.seq) and int(predicted_region[i].split()[0]) != 0 : #which means it already reads target lines and binds to target sRNA at that region
            predicted_base_pairs.append(int(predicted_region[i].split()[1]))
        else:
            pass

    predicted_base_pairs=set(predicted_base_pairs)

    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))

    if True_Positives != 0 or False_Positives != 0:
        print_results(args.program,True_Positives,False_Positives,False_Negatives)
    else:
        print_results(args.program,True_Positives,1,False_Negatives) #add pseudo count for bifold

    os.remove(RNA_temp_file.name)
    os.remove(sRNA_temp_file.name)

    return


def DuplexFold_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()


    sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='DuplexFold_sRNA_',suffix='.fasta',mode="w+",delete=False)
    sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.seq)))
    sRNA_temp_file.close()

    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='DuplexFold_RNA_',suffix='.fasta',mode="w+",delete=False)
    RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq)))
    RNA_temp_file.close()

    shell_command="""DuplexFold -m 1 %s %s %s > /dev/null && cat %s | grep -v ENERGY | gawk '{print $(NF-1),$NF}'""" %(sRNA_temp_file.name,RNA_temp_file.name,RNA_temp_file.name+".out",RNA_temp_file.name+".out")

    predicted_region=subprocess.check_output(shell_command,shell=True).splitlines()

    true_base_pairs=set(range(real_region[0],real_region[1]+1))

    predicted_base_pairs=list()
    for i in range(0,len(predicted_region)):
        if i+1 > len(sRNA.seq) and int(predicted_region[i].split()[0]) < len(sRNA.seq) and int(predicted_region[i].split()[0]) != 0: #which means it already reads target lines and binds to target sRNA at that region
            predicted_base_pairs.append(int(predicted_region[i].split()[1]))
        else:
            pass

    predicted_base_pairs=set(predicted_base_pairs)

    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))

    print_results(args.program,True_Positives,False_Positives,False_Negatives)

    os.remove(RNA_temp_file.name)
    os.remove(sRNA_temp_file.name)

    return

'''
# SSEARCH 36.3.6 Jan, 2014(preload9)
# Query: HBII-239 - 86 nt
# Database: 5.8S.human.fasta
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 2 hits found
HBII-239	gb-U13369.1	66.67	48	14	2	40	86	2	48	0.033	18.6
HBII-239	gb-U13369.1	52.44	82	33	6	81	1	36	112	1	11.0
# SSEARCH processed 1 queries
'''
def ssearch_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()

    sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='ssearch_sRNA_',suffix='.fasta',mode="w+",delete=False)
    sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.reverse_complement().seq)))
    sRNA_temp_file.close()

    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='ssearch_RNA_',suffix='.fasta',mode="w+",delete=False)
    RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq)))
    RNA_temp_file.close()


    shell_command="""ssearch36 -s '/home/suu13/rna.mat' %s %s -m 8C | grep -v '#' """ %(sRNA_temp_file.name,RNA_temp_file.name)

    predicted_region=subprocess.check_output(shell_command,shell=True).splitlines() #returns only start of target RNA

    predicted_base_pairs=list()

    for i in predicted_region:
        predicted_base_pairs=predicted_base_pairs+ range(int(i.split()[8]),int(i.split()[9])+1)

    predicted_base_pairs=set(predicted_base_pairs)

    true_base_pairs=set(range(real_region[0],real_region[1]+1))

    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))


    print_results(args.program,True_Positives,False_Positives,False_Negatives)

    os.remove(RNA_temp_file.name)
    os.remove(sRNA_temp_file.name)
    return


def ractip_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()

    sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='ractip_sRNA_',suffix='.fasta',mode="w+",delete=False)
    sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.reverse_complement().seq)))
    sRNA_temp_file.close()

    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='ractip_RNA_',suffix='.fasta',mode="w+",delete=False)
    RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq)))
    RNA_temp_file.close()

    shell_command="""ractip %s %s | tail -n 1""" %(sRNA_temp_file.name,RNA_temp_file.name)


    predicted_region=subprocess.check_output(shell_command,shell=True)


    predicted_base_pairs=set(map(lambda  x: x+1,[i.start() for i in finditer(']',predicted_region)])) # add +1 too all index locations and convert to a set

    true_base_pairs=set(range(real_region[0],real_region[1]+1))

    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))


    print_results(args.program,True_Positives,False_Positives,False_Negatives)

    os.remove(RNA_temp_file.name)
    os.remove(sRNA_temp_file.name)
    return


def bistarna_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()

    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='bistarna_RNA_',suffix='.fasta',mode="w+",delete=False)
    RNA_temp_file.write(">%s\n%s\n>%s\n%s\n" %(str(sRNA.description),str(sRNA.seq),str(RNA.description),str(RNA.seq)))
    RNA_temp_file.close()

    shell_command="""bistarna %s | grep ']'""" %(RNA_temp_file.name)


    predicted_region=subprocess.check_output(shell_command,shell=True)


    predicted_base_pairs=set(map(lambda  x: x+1,[i.start() for i in finditer(']',predicted_region)])) # add +1 too all index locations of predicted regions and convert to a set

    true_base_pairs=set(range(real_region[0],real_region[1]+1))

    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))


    print_results(args.program,True_Positives,False_Positives,False_Negatives)

    os.remove(RNA_temp_file.name)

    return

def AccessFold_benchmark(Input_sRNA_File,Input_Target_File,real_region):
    sRNA=SeqIO.parse(Input_sRNA_File,"fasta").next()
    RNA=SeqIO.parse(Input_Target_File,"fasta").next()


    sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='AccessFold_sRNA_',suffix='.fasta',mode="w+",delete=False)
    sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.seq)))
    sRNA_temp_file.close()

    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='AccessFold_RNA_',suffix='.fasta',mode="w+",delete=False)
    RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq)))
    RNA_temp_file.close()

    shell_command="""AccessFold -m 1 %s %s %s > /dev/null && cat %s | grep -v ENERGY | gawk '{print $(NF-1),$NF}'""" %(sRNA_temp_file.name,RNA_temp_file.name,RNA_temp_file.name+".out",RNA_temp_file.name+".out")


    predicted_region=subprocess.check_output(shell_command,shell=True).splitlines()

    true_base_pairs=set(range(real_region[0],real_region[1]+1))

    predicted_base_pairs=list()
    for i in range(0,len(predicted_region)):
        if i+1 > len(sRNA.seq) and int(predicted_region[i].split()[0]) < len(sRNA.seq) and int(predicted_region[i].split()[0]) != 0: #which means it already reads target lines and binds to target sRNA at that region
            predicted_base_pairs.append(int(predicted_region[i].split()[1]))
        else:
            pass

    predicted_base_pairs=set(predicted_base_pairs)



    True_Positives=len(true_base_pairs.intersection(predicted_base_pairs))
    False_Positives=len(predicted_base_pairs.difference(true_base_pairs))
    False_Negatives=len(true_base_pairs.difference(predicted_base_pairs))

    if True_Positives != 0 or False_Positives != 0:
        print_results(args.program,True_Positives,False_Positives,False_Negatives)
    else:
        print_results(args.program,True_Positives,1,False_Negatives) #add pseudo count for AccessFold


    os.remove(RNA_temp_file.name)
    os.remove(sRNA_temp_file.name)

    return



def Parse_Dot_Bracket(secondary_structure):
    stack=list()
    predicted_base_pairs=list()
    for i in range(0,len(secondary_structure)):
        if len(stack) > 0:
            if stack[-1] == '(' and secondary_structure[i] == ')':
                stack.pop()
            elif stack[-1] == '(' and secondary_structure[i] == '(':
                stack.append(secondary_structure[i])
            elif secondary_structure[i]  == ')':
                predicted_base_pairs.append(i+1)
            elif secondary_structure[i] == '.':
                pass
        elif secondary_structure[i] == '(':
            stack.append(secondary_structure[i])
        elif secondary_structure[i] == ')':
            predicted_base_pairs.append(i+1)
        else:
            pass
    return predicted_base_pairs


def print_results(Algorithm,True_Positives,False_Positives,False_Negatives):

    TPR=float(True_Positives/float(True_Positives+False_Negatives))
    PPV=float(True_Positives/float(True_Positives+False_Positives))
    MCC=sqrt(TPR*PPV)
    #F=float(2*PPV*TPR/float(PPV+TPR))

    print "Algorithm\tTrue_Positives\tFalse_Positives\tFalse_Negatives\tTPR\tPPV\tMCC"
    print "%s\t%d\t%d\t%d\t%f\t%f\t%f\t" % (Algorithm,True_Positives,
                                      False_Positives,False_Negatives,
                                      TPR,PPV,MCC)

    return


def main():

    #home-made switch
    program_function={'RIsearch': lambda: RIsearch_benchmark(args.sRNA,args.targetRNA,args.window),
                      'IntaRNA': lambda: IntaRNA_benchmark(args.sRNA,args.targetRNA,args.window),
                      'RNAplex': lambda: RNAplex_benchmark(args.sRNA,args.targetRNA,args.window),
                      'RNAcofold': lambda: RNAcofold_benchmark(args.sRNA,args.targetRNA,args.window),
                      'pairfold': lambda: Pairfold_benchmark(args.sRNA,args.targetRNA,args.window),
                      'RNAup': lambda: RNAup_benchmark(args.sRNA,args.targetRNA,args.window),
                      'RNAduplex': lambda: RNAduplex_benchmark(args.sRNA,args.targetRNA,args.window),
                      'RNAhybrid': lambda: RNAhybrid_benchmark(args.sRNA,args.targetRNA,args.window),
                      'bifold': lambda: bifold_benchmark(args.sRNA,args.targetRNA,args.window),
                      'DuplexFold': lambda: DuplexFold_benchmark(args.sRNA,args.targetRNA,args.window),
                      'ssearch': lambda: ssearch_benchmark(args.sRNA,args.targetRNA,args.window),
                      'ractip': lambda: ractip_benchmark(args.sRNA,args.targetRNA,args.window),
                      'bistarna': lambda: bistarna_benchmark(args.sRNA,args.targetRNA,args.window),
                      'AccessFold': lambda: AccessFold_benchmark(args.sRNA,args.targetRNA,args.window)

                        }


    program_function[args.program]()






if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="RNA_interaction_wrapper.py")
    Argument_Parser.add_argument('-program',type=str,help="Program to run",choices=['RIsearch','IntaRNA','RNAplex','RNAcofold',
                                                                                    'pairfold','RNAup','RNAduplex','RNAhybrid','bifold','DuplexFold','ssearch','ractip',
                                                                                    'bistarna','AccessFold'],required=True)
    Argument_Parser.add_argument('-sRNA',type=str,help="Small RNA file/miRNA etc.",required=True)
    Argument_Parser.add_argument('-targetRNA',type=str,help="Target RNA file/mRNAs/UTRs etc.",required=True)
    Argument_Parser.add_argument('-window',type=int,nargs=2,help="Target interaction region, from start to stop",required=True)
    args=Argument_Parser.parse_args()
    main()

