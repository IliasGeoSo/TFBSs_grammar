from scipy.stats import binom_test
from random import shuffle
import numpy as np

def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

    

counted=0
same_sig=0;opposite_sig=0;
files=glob.glob("inter_prom/*")
datafile=open("Asym_TF_same_opposite_genomic_prom.txt","w")
datafile.write("TF"+'\t'+"strand asymmetry"+'\t'+"p_val"+'\n')
ScoresL=[]
ScoresL_std=[]
for one in files:
        print(one)
        DataL=reader(one)
        TF=DataL[0][9]
        same=0;opposite=0;
        if len(DataL)>1000:
                for k in range(0,len(DataL)-1):
                        end1=int(DataL[k][2])
                        start2=int(DataL[k+1][1])
                        if start2-end1<50 and end1<start2:
                                if DataL[k][-3]==DataL[k+1][-3]:
                                        if DataL[k][3]==DataL[k+1][3]:
                                                same+=1
                                        if DataL[k][3]!=DataL[k+1][3] and DataL[k][3] in ["+","-"] and DataL[k+1][3] in ["+","-"]:
                                                opposite+=1

                Same_simL=[];Opposite_simL=[];
                for sim in range(10):
                        Strands=[mn[-3] for mn in DataL]
                        shuffle(Strands)
                        same_sim=0;opposite_sim=0
                        for k in range(0,len(DataL)-1):
                                        end1=int(DataL[k][2])
                                        start2=int(DataL[k+1][1])
                                        if start2-end1<50 and end1<start2:
                                                if Strands[k]==Strands[k+1]:
                                                        same_sim+=1
                                                if Strands[k]!=Strands[k+1] and Strands[k] in ["+","-"] and Strands[k+1] in ["+","-"]:
                                                        opposite_sim+=1
                        Same_simL+=[same_sim];Opposite_simL+=[opposite_sim];
                if np.mean(Same_simL)+np.mean(Opposite_simL)!=0 and float(np.mean(Same_simL))/float(np.mean(Same_simL)+np.mean(Opposite_simL))!=0 and same+opposite!=0:
                        ScoresL+=[(float(same/float(same+opposite))/float(np.mean(Same_simL)))/float(np.mean(Same_simL)+np.mean(Opposite_simL))]

        p_val=min(1,binom_test(same,same+opposite,0.5)*845)
        if same+opposite!=0:
                datafile.write(TF+'\t'+str(same/float(same+opposite))+'\t'+str(p_val)+'\n')
        if binom_test(same,same+opposite,0.5)*835<0.05:
                counted+=1
                if same>opposite:
                        same_sig+=1

print(counted,same_sig,len(files))
ScoresL, list2 = zip(*sorted(zip(list1, list2)))
print(counted,same_sig,len(files))
