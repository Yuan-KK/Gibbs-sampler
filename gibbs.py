#!/usr/bin/python3
# -*- coding: UTF-8 -*-

################################################
# Finding motif using Gibbs sampling by Keke Yuan (1452915529@qq.com)
# Date: Sep. 25, 2021  
################################################

import datetime, argparse
from os.path import exists
from os import mkdir
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import seqlogo
from matplotlib.backends.backend_pdf import PdfPages

def get_opt():
	'''Check and parsing the opts'''
	parser=argparse.ArgumentParser(usage='python gibbs.py -i read.fasta -L (int) [-o output/]')
	parser.add_argument('-i','--input',nargs='?',type=str,default=sys.stdin, help="The path of fasta file, the suffix must be '.fa' or '.fasta'. [Required]",required=True)
	parser.add_argument('-o','--output',default = './output', help='output folder[Required]',required=True)
	parser.add_argument('-L','--length',nargs='?',type=int,default='query',help="The length of motif expected to get",required=True)
	args=parser.parse_args()
	return args

def seqcut(df,startlist,L):
    seqcutlist = []
    for i in range(len(startlist)):
        s = startlist[i]
        h = df['seq'][i]
        d = list(h)[s:s+L]
        seqcutlist.append(d)
    seqcut_df = pd.DataFrame(seqcutlist)
    return seqcut_df
    
# Count matrix
def cmatrix(seqcut_df,L,nA,nC,nG,nT):
    i = 0 
    dicc = {}
    while i < L:
        dicc.setdefault('A',[]).append(list(seqcut_df[i]).count('A')) 
        dicc.setdefault('C',[]).append(list(seqcut_df[i]).count('C'))
        dicc.setdefault('G',[]).append(list(seqcut_df[i]).count('G'))
        dicc.setdefault('T',[]).append(list(seqcut_df[i]).count('T'))
        i = i+1

    Fa = nA - sum(dicc['A'])
    Fc = nC - sum(dicc['C'])
    Fg = nG - sum(dicc['G'])
    Ft = nT - sum(dicc['T'])

    dicc["A"].insert(0,Fa)
    dicc["C"].insert(0,Fc)
    dicc["G"].insert(0,Fg)
    dicc["T"].insert(0,Ft)
    return dicc

# PPM matrix
def PPM(matrix,L,a=0.0001):  # Cmatrix
    basec = matrix["A"][0]+matrix["C"][0]+matrix["G"][0]+matrix["T"][0]
    nc = matrix["A"][1]+matrix["C"][1]+matrix["G"][1]+matrix["T"][1]
    i = 1
    PPM = {}
    while i <= L:
        PPM.setdefault('A',[]).append((matrix["A"][i]+a*matrix["A"][0]/basec)*basec / (nc*matrix["A"][0]))
        PPM.setdefault('C',[]).append((matrix["C"][i]+a*matrix["C"][0]/basec)*basec / (nc*matrix["C"][0]))
        PPM.setdefault('G',[]).append((matrix["G"][i]+a*matrix["G"][0]/basec)*basec / (nc*matrix["G"][0]))
        PPM.setdefault('T',[]).append((matrix["T"][i]+a*matrix["T"][0]/basec)*basec / (nc*matrix["T"][0]))
        i = i+1
    return PPM
    

# F value
def F(C,P,L):
    C = pd.DataFrame.from_dict(C,orient='index')
    P = pd.DataFrame.from_dict(P,orient='index')
    Fvalue = 0
    for i in range(4):
        for j in range(L):
            Fvalue += C.iat[i,j+1]*P.iat[i,j]
    return Fvalue

# Score matrix
def cutsample(df,N,L):  
    i = 0 
    countlist = []
    while i < L:
        pera = (list(df[i]).count('A')+1) / (N + 4)   
        perc = (list(df[i]).count('C')+1) / (N + 4)
        perg = (list(df[i]).count('G')+1) / (N + 4)
        pert = (list(df[i]).count('T')+1) / (N + 4)
        countlist.append([pera,perc,perg,pert])
        i = i+1
        if i == L:
            break
    countlist = np.asarray(countlist).T
    return countlist

# Print the Maximum of one list.
# If there are multiple Maximums, select one of them at random.
def mutimaxindex(nums):
    max_of_nums = max(nums)
    tup = [(i, nums[i]) for i in range(len(nums))]
    maxlist = [i for i, n in tup if n == max_of_nums]
    a = random.choice(maxlist)
    return(a)

# Outputs the number of occurrences of a repeated value.
def continue_num(lst):
    length = len(lst)
    total_num = []
    j = 1
    for i in range(length - 1):
        if lst[i] == lst[i+1]:
            j += 1
        else:
            total_num.append(j)
            j = 1
    total_num.append(j)
    fremax = max(total_num)
    return(fremax)
    
# One loop
def gibbs(dna_df,N,L,nA,nC,nG,nT,t=1000):
    i = 0
    # Random selection of starting sites
    startlist = []
    for i in range(N):
        l = len(dna_df['seq'][i])
        s = random.randint(0,l-L)  
        i += 1
        startlist.append(s)
    # Initial motif matrix
    motifdf = seqcut(dna_df,startlist,L)
    c0matrix = cmatrix(motifdf,L,nA,nC,nG,nT)
    p0matrix = PPM(c0matrix,L)
    # Initial F value
    F1 = F(c0matrix,p0matrix,L)   
    Flist = []
    Flist.append(F1)
    u = 1
    df_o = motifdf
    i = 0
    h = [x for x in range(N)]
    random.shuffle(h)
    while 1:   
        for i in h:
            df_S = dna_df.loc[i]
            df_S = list(df_S[0]) 
            df_new = motifdf.drop([i], axis=0)
            xmatrix = cutsample(df_new,N-1,L)
            kmers = []
            m = 0 
            while m <= len(df_S)-L:
                n = m + L
                di = df_S[m:n]
                di = ''.join(di)
                kmers.append(di)
                m += 1 
                if m > len(df_S)-L:
                    break
            x = 0 
            y = 0
            scorelist = []
            k = len(kmers)
            for x in range(k) :
                oddratio = 1
                for y in range(1,L):
                    if list(kmers[x])[y] == 'A':   
                        odd = xmatrix[0,y]
                    elif list(kmers[x])[y] == 'C':
                        odd = xmatrix[1,y]
                    elif list(kmers[x])[y] == 'G':
                        odd = xmatrix[2,y]
                    elif list(kmers[x])[y] == 'T':
                        odd = xmatrix[3,y]
                    y += 1
                    oddratio *= odd
                scorelist.append(oddratio)
                x += 1
            maxindex = mutimaxindex(scorelist)
            startlist[i] = maxindex
            motifdf.loc[i] = list(dna_df['seq'][i])[maxindex:maxindex+L] 
            c1matrix = cmatrix(motifdf,L,nA,nC,nG,nT)
            p1matrix = PPM(c1matrix,L)
            F2 = F(c1matrix,p1matrix,L)
            Flist.append(F2) 
            u += 1    
            if i == h[-1]:
                h = [x for x in range(N)]
                random.shuffle(h)
                if F2 >=  F1:
                    F1 = F2
                    tstartlist = {}
                    tseqcut={}
                    tF={}
                    top = [x for x in range(-L,L+1)]
                    for j in top:
                        b = [x+j for x in startlist]
                        tstartlist[j]=b
                        tdf=seqcut(dna_df,b,L) 
                        tseqcut[j] = tdf
                        tcmatrix = cmatrix(tdf,L,nA,nC,nG,nT)
                        tpmatrix = PPM(tcmatrix,L)
                        tF[j]= F(tcmatrix,tpmatrix,L)
                    tmaxindex = max(tF, key=tF.get)
                    startlist = tstartlist[tmaxindex]
                    motifdf = tseqcut[tmaxindex]
                    df_o= motifdf.copy(deep=True)
                    F2 = tF[tmaxindex]
                else:
                    motifdf = df_o.copy(deep=True)
        if u > t:
            break
    return F2,motifdf,Flist

def main(args):
    input = args.input
    output = args.output.rstrip('/')
    if not exists(output):
        mkdir(output)
    L = args.length

    log_file=open(output+'/outfile.log',"a") 
    log_file.write("Start: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M") + "\n")  
    
    sys.stdout=log_file
    pd.set_option('display.max_columns', None) 
    pd.set_option('display.max_rows', None)

    with open(input,'r') as f1:
        with open(output+'/tmp.txt','w') as f2:
            for line in f1:
                if 'A' in line:
                    line=line.strip('\n')
                f2.writelines(line)
    linelist = []
    with open(output+'/tmp.txt','r') as f2:
        for line in f2:
            if '>' not in line:
                linelist.append(line)
    for n in range(len(linelist)) :
        if linelist[n].endswith('\n'):
            linelist[n]=linelist[n][:-1]
    dna_seq = {'seq':linelist}
    metadna_df = pd.DataFrame(dna_seq)  
    N = len(linelist)    
    nA = 0
    nC = 0
    nG = 0
    nT = 0
    for i in range(N):
        nA += metadna_df['seq'][i].count("A")
        nC += metadna_df['seq'][i].count("C")
        nG += metadna_df['seq'][i].count("G")
        nT += metadna_df['seq'][i].count("T")

    diciri_F={}
    diciri_df={}
    F2list = []
    ircount = 0
    with PdfPages(output+"/plot.pdf") as pdf:
        while 1:
            ircount += 1
            diciri_F[ircount],diciri_df[ircount],Flist=gibbs(metadna_df,N,L,nA,nC,nG,nT,t=1000)
            plt.title(ircount, fontsize = 25)
            plt.plot(Flist)
            pdf.savefig()
            plt.close()
            imaxindex = max(diciri_F, key=diciri_F.get)
            F2_max=diciri_F[imaxindex]
            if diciri_F[ircount] >= F2_max :
                F2list.append(diciri_F[ircount])
                max_num = continue_num(F2list)
                if max_num >= 2:
                    break

    print("Gibbs sampler reached convergence after "+str(ircount)+" iterations.","\n")
    imaxindex = max(diciri_F, key=diciri_F.get)
    F2=diciri_F[imaxindex]
    motif_df= diciri_df[imaxindex]
    ff = pd.DataFrame(cmatrix(motif_df,L,nA,nC,nG,nT))
    pfm = ff.drop([0], axis=0)
    cpm = seqlogo.CompletePm(pfm)
    print("Global alignment score (F) = ",F2,"\n")
    print("Motif"+"\n",motif_df,"\n")
    print("Final site-specific counts:\n",pfm,"\n")
    print("Final PWM:\n",cpm.pwm,"\n")
    with open(output+"/motif.csv",'w',newline='') as f1:
        motif_df.to_csv(f1,index=F)

    with open(output+"/matrix.csv",'w',newline='') as f2:
        matrixdf=cutsample(motif_df,N,L)
        matrixdf = pd.DataFrame(matrixdf,index = ['A','C','G','T'])
        matrixdf.to_csv(f2)

    log_file.write("Finish: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M") + '\n')
    log_file.close()

if __name__=='__main__':
	args=get_opt()
	main(args)
