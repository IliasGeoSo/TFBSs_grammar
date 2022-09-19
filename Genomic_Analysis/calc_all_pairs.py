import re,os,sys,glob,math
from collections import Counter
from scipy.stats import binom_test

def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[] 
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

Dic={}
Dic_strands={}
DataL=reader("All_motifs.bed")
for k in range(0,len(DataL)-1,1):
        end1=int(DataL[k][2])
        start2=int(DataL[k+1][2])
        counter=0
        while end1<start2 and abs(end1-start2)<100 and k+counter+1<len(DataL):
                counter+=1
                TF1=DataL[k][-1]
                TF2=DataL[k+1][-1]
                try:
                        Dic[TF1]+=[TF2]
                except:
                        Dic[TF1]=[TF2]

                try:
                        Dic[TF2]+=[TF1]
                except:
                        Dic[TF2]=[TF1]

                try:
                        Dic_strands[TF1+"-"+TF2]+=DataL[k][3]+DataL[k+1][3]
                except:
                        Dic_strands[TF1+"-"+TF2]=DataL[k][3]+DataL[k+1][3]

                try:
                        Dic_strands[TF2+"-"+TF1]+=DataL[k+1][3]+DataL[k][3]
                except:
                        Dic_strands[TF2+"-"+TF1]=DataL[k+1][3]+DataL[k][3]
                if k+counter+1<len(DataL):
                        start2=int(DataL[k+1+counter][2])
        else:
                #end1=int(DataL[k][2])
                #start2=int(DataL[k+1][2])
                counter=0
                pass

biased=0;biased2=0;
for TF in Dic.keys():
        Total1=len(Dic[TF])
        occs=Counter(Dic[TF])
        TFspairs=[];jacsL=[];
        for k in occs.items():
                jac=k[1]/float(max(Total1,len(Dic[k[0]])))
                jacsL+=[jac]
                TFspairs+=[TF+"-"+k[0]]

        jacsL, TFspairs = zip(*sorted(zip(jacsL, TFspairs)))
        #print(jacsL[-10:],TFspairs[-10:])
        jacsL=jacsL[-10:];TFspairs=TFspairs[-10:]
        for TF_pair_used in TFspairs:
                same=Dic_strands[TF_pair_used].count("++")+Dic_strands[TF_pair_used].count("--")
                opposite=Dic_strands[TF_pair_used].count("+-")+Dic_strands[TF_pair_used].count("-+")
                if binom_test(same,same+opposite,0.5)*len(Dic_strands.keys())<0.05:
                        biased+=1

print(biased/float(len(Dic_strands.keys())))

