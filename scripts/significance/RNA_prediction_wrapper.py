#!/usr/bin/env python2.7


'''
Created on 27/10/2015

@author: suu13
'''


import argparse
import re
from Bio import SeqIO
import subprocess
import tempfile
import os
from joblib import Parallel,delayed
import time


__author__ = 'suu13'



def RNAup_Execute_Parallel(Input_sRNA_File,Input_Target_File,real_region):
    sRNA_FASTA=list(SeqIO.parse(Input_sRNA_File,"fasta"))
    Target_FASTA=list(SeqIO.parse(Input_Target_File,"fasta"))

    if real_region is not None:
        size=abs(real_region[1]-real_region[0])+5 # for bacteria +5 is enough
    else:
        size=25 #RNAup default


    for sRNA in sRNA_FASTA:
        sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='RNAup_sRNA_',suffix='.tmp',mode="w+",delete=False)
        sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.seq)))
        sRNA_temp_file.close()

        RNAjobList=[]

        while True:
        #for RNA in Target_FASTA:
            try:
                if len(RNAjobList) < args.cpu:
                    RNA=Target_FASTA.pop(0)

                    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='RNAup_RNA_',suffix='.tmp',mode="w+",delete=False)
                    RNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.seq))) #write sRNA to the first line
                    RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq))) #write RNA to the second line
                    RNA_temp_file.close()
                    shell_command="""~/progs/ViennaRNA-2.1.9/Progs/RNAup -w %d -b -o --interaction_first < %s | gawk 'match($0,/ \((.*) =/,m) {printf m[1]}'""" %(size,RNA_temp_file.name)
                    job_id=subprocess.Popen(shell_command,shell=True,stdout=subprocess.PIPE)
                    RNAjobList.append([str(sRNA.description),str(RNA.description),job_id,RNA_temp_file.name])
                else:
                    time.sleep(0.1) #wait for 0.1s to check status of jobs
                    for job in RNAjobList:
                        if subprocess.Popen.poll(job[2])!=None: #check the status of job object
                            os.remove(job[3])
                            print "%s\t%s\t%s" % (job[0],job[1],job[2].communicate()[0].strip())
                            RNAjobList.remove(job)
            except:
                    while len(RNAjobList) != 0:
                        time.sleep(0.1) #wait for 0.1s to check status of jobs
                        for job in RNAjobList:
                            if subprocess.Popen.poll(job[2])!=None: #check the status of job object
                                os.remove(job[3])
                                print "%s\t%s\t%s" % (job[0],job[1],job[2].communicate()[0].strip())
                                RNAjobList.remove(job)
                    break


        os.remove(sRNA_temp_file.name)

    return


def Pairfold_Execute(Input_sRNA_File,Input_Target_File):
    sRNA_FASTA=list(SeqIO.parse(Input_sRNA_File,"fasta"))
    Target_FASTA=list(SeqIO.parse(Input_Target_File,"fasta"))

    for sRNA in sRNA_FASTA:
        for RNA in Target_FASTA:
            shell_command="""/home/suu13/misc_stuff/MultiRNAFold-2.0/pairfold "%s" "%s" -m RNA | grep MFE | awk '{print $NF}' """ % (str(sRNA.seq),str(RNA.seq)) #enter pairfold actual path
            result=subprocess.check_output(shell_command,shell=True)
            print "%s\t%s\t%s" % (str(sRNA.description),str(RNA.description),result.strip())
    return



def Pairfold_Execute_Parallel(Input_sRNA_File,Input_Target_File):
    sRNA_FASTA=list(SeqIO.parse(Input_sRNA_File,"fasta"))
    Target_FASTA=list(SeqIO.parse(Input_Target_File,"fasta"))

    for sRNA in sRNA_FASTA:
        RNAjobList=[]

        #for RNA in Target_FASTA:
        while True:
            try:
                if len(RNAjobList) < args.cpu:
                    RNA=Target_FASTA.pop(0)

                    shell_command="""/home/suu13/misc_stuff/MultiRNAFold-2.0/pairfold "%s" "%s" -m RNA | grep MFE | awk '{print $NF}' """ % (str(sRNA.seq),str(RNA.seq)) #enter pairfold actual path
                    job_id=subprocess.Popen(shell_command,shell=True,stdout=subprocess.PIPE)
                    RNAjobList.append([str(sRNA.description),str(RNA.description),job_id])
                else:
                    time.sleep(0.1) #wait for 0.1s to check status of jobs
                    for job in RNAjobList:
                        if subprocess.Popen.poll(job[2])!=None: #check the status of job object
                            print "%s\t%s\t%s" % (job[0],job[1],job[2].communicate()[0].strip())
                            RNAjobList.remove(job)
            except:
                    while len(RNAjobList) != 0:
                        time.sleep(0.1) #wait for 0.1s to check status of jobs
                        for job in RNAjobList:
                            if subprocess.Popen.poll(job[2])!=None: #check the status of job object
                                print "%s\t%s\t%s" % (job[0],job[1],job[2].communicate()[0].strip())
                                RNAjobList.remove(job)
                    break



    return


def bifold_Execute(Input_sRNA_File,Input_Target_File):
    sRNA_FASTA=list(SeqIO.parse(Input_sRNA_File,"fasta"))
    Target_FASTA=list(SeqIO.parse(Input_Target_File,"fasta"))

    for sRNA in sRNA_FASTA:
        sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='bifold_sRNA_',suffix='.tmp',mode="w+",delete=False)
        sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.seq)))
        sRNA_temp_file.close()
        for RNA in Target_FASTA:
            RNA_temp_file=tempfile.NamedTemporaryFile(prefix='bifold_RNA_',suffix='.tmp',mode="w+",delete=False)
            RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq)))
            RNA_temp_file.close()
            shell_command="""bifold -m 1 %s %s %s > /dev/null && head -n 1 %s | gawk '{match($0,/ENERGY = (.*) /,m); print m[1]}' """ %(sRNA_temp_file.name,RNA_temp_file.name,RNA_temp_file.name+".out",RNA_temp_file.name+".out")
            result=subprocess.check_output(shell_command,shell=True)
            print "%s\t%s\t%s" % (str(sRNA.description),str(RNA.description),result.strip())
            os.remove(RNA_temp_file.name)
        os.remove(sRNA_temp_file.name)

    return


def bifold_Execute_Parallel(Input_sRNA_File,Input_Target_File):
    sRNA_FASTA=list(SeqIO.parse(Input_sRNA_File,"fasta"))
    Target_FASTA=list(SeqIO.parse(Input_Target_File,"fasta"))

    for sRNA in sRNA_FASTA:
        sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='bifold_sRNA_',suffix='.tmp',mode="w+",delete=False)
        sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.seq)))
        sRNA_temp_file.close()

        RNAjobList=[]

        while True:
        #for RNA in Target_FASTA:
            try:
                if len(RNAjobList) < args.cpu:
                    RNA=Target_FASTA.pop(0)

                    RNA_temp_file=tempfile.NamedTemporaryFile(prefix='bifold_RNA_',suffix='.tmp',mode="w+",delete=False)
                    RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq)))
                    RNA_temp_file.close()
                    shell_command="""bifold -m 1 %s %s %s > /dev/null && head -n 1 %s | gawk '{match($0,/ENERGY = (.*) /,m); print m[1]}' """ %(sRNA_temp_file.name,RNA_temp_file.name,RNA_temp_file.name+".out",RNA_temp_file.name+".out")
                    #result=subprocess.check_output(shell_command,shell=True)
                    job_id=subprocess.Popen(shell_command,shell=True,stdout=subprocess.PIPE)
                    #print "%s\t%s\t%s" % (str(sRNA.description),str(RNA.description),result.strip())
                    RNAjobList.append([str(sRNA.description),str(RNA.description),job_id,RNA_temp_file.name])
                    #os.remove(RNA_temp_file.name)
                else:
                    time.sleep(0.1) #wait for 0.1s to check status of jobs
                    for job in RNAjobList:
                        if subprocess.Popen.poll(job[2])!=None: #check the status of job object
                            os.remove(job[3])
                            print "%s\t%s\t%s" % (job[0],job[1],job[2].communicate()[0].strip())
                            RNAjobList.remove(job)
            except:
                    while len(RNAjobList) != 0:
                        time.sleep(0.1) #wait for 0.1s to check status of jobs
                        for job in RNAjobList:
                            if subprocess.Popen.poll(job[2])!=None: #check the status of job object
                                os.remove(job[3])
                                print "%s\t%s\t%s" % (job[0],job[1],job[2].communicate()[0].strip())
                                RNAjobList.remove(job)
                    break


        os.remove(sRNA_temp_file.name)

    return

def DuplexFold_Execute(Input_sRNA_File,Input_Target_File):
    sRNA_FASTA=list(SeqIO.parse(Input_sRNA_File,"fasta"))
    Target_FASTA=list(SeqIO.parse(Input_Target_File,"fasta"))

    for sRNA in sRNA_FASTA:
        sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='DuplexFold_sRNA_',suffix='.tmp',mode="w+",delete=False)
        sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.seq)))
        sRNA_temp_file.close()
        for RNA in Target_FASTA:
            RNA_temp_file=tempfile.NamedTemporaryFile(prefix='DuplexFold_RNA_',suffix='.tmp',mode="w+",delete=False)
            RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq)))
            RNA_temp_file.close()
            shell_command="""DuplexFold -m 1 %s %s %s > /dev/null && head -n 1 %s | gawk '{match($0,/ENERGY = (.*) /,m); print m[1]}' """ %(sRNA_temp_file.name,RNA_temp_file.name,RNA_temp_file.name+".out",RNA_temp_file.name+".out")
            result=subprocess.check_output(shell_command,shell=True)
            print "%s\t%s\t%s" % (str(sRNA.description),str(RNA.description),result.strip())
            os.remove(RNA_temp_file.name)
        os.remove(sRNA_temp_file.name)

    return


def AccessFold_Execute(Input_sRNA_File,Input_Target_File):
    sRNA_FASTA=list(SeqIO.parse(Input_sRNA_File,"fasta"))
    Target_FASTA=list(SeqIO.parse(Input_Target_File,"fasta"))

    for sRNA in sRNA_FASTA:
        sRNA_temp_file=tempfile.NamedTemporaryFile(prefix='AccessFold_sRNA_',suffix='.tmp',mode="w+",delete=False)
        sRNA_temp_file.write(">%s\n%s\n" %(str(sRNA.description),str(sRNA.seq)))
        sRNA_temp_file.close()
        for RNA in Target_FASTA:
            RNA_temp_file=tempfile.NamedTemporaryFile(prefix='AccessFold_RNA_',suffix='.tmp',mode="w+",delete=False)
            RNA_temp_file.write(">%s\n%s\n" %(str(RNA.description),str(RNA.seq)))
            RNA_temp_file.close()
            shell_command="""AccessFold -m 1 %s %s %s > /dev/null && head -n 1 %s | gawk '{match($0,/ENERGY = (.*) /,m); print m[1]}' """ %(sRNA_temp_file.name,RNA_temp_file.name,RNA_temp_file.name+".out",RNA_temp_file.name+".out")
            result=subprocess.check_output(shell_command,shell=True)
            print "%s\t%s\t%s" % (str(sRNA.description),str(RNA.description),result.strip())
            os.remove(RNA_temp_file.name)
        os.remove(sRNA_temp_file.name)

    return




def main():

    #home-made switch
    program_function={ 'pairfold': lambda: Pairfold_Execute(args.sRNA,args.targetRNA) if args.cpu is None else Pairfold_Execute_Parallel(args.sRNA,args.targetRNA),
                       'bifold': lambda: bifold_Execute(args.sRNA,args.targetRNA) if args.cpu is None else bifold_Execute_Parallel(args.sRNA,args.targetRNA),
                       'DuplexFold': lambda: DuplexFold_Execute(args.sRNA,args.targetRNA),
                       'RNAup': lambda: RNAup_Execute_Parallel(args.sRNA,args.targetRNA,args.window),
                       'AccessFold': lambda: AccessFold_Execute(args.sRNA,args.targetRNA),
                        }


    program_function[args.program]()







if __name__ == '__main__':
    Argument_Parser=argparse.ArgumentParser(prog="RNA_prediction_wrapper.py")
    Argument_Parser.add_argument('-program',type=str,help="Program to run",choices=['pairfold','bifold','DuplexFold','RNAup','AccessFold'],required=True)
    Argument_Parser.add_argument('-sRNA',type=str,help="Small RNA file/miRNA etc.",required=True)
    Argument_Parser.add_argument('-targetRNA',type=str,help="Target RNA file/mRNAs etc.",required=True)
    Argument_Parser.add_argument('-cpu',type=int,help="Use parallel computing")
    Argument_Parser.add_argument('-window',type=int,nargs=2,help="Target interaction region, from start to stop",required=False)
    args=Argument_Parser.parse_args()
    main()
