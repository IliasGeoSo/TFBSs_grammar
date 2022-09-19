import re,os,sys,glob,random,pdb,math
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from scipy.stats import pearsonr,ttest_ind,kruskal,spearmanr,f_oneway
from scipy.stats import wilcoxon
import numpy as np
from Bio.Seq import Seq
from matplotlib.pyplot import figure
import pandas as pd
plt.rcParams.update({'font.size': 16})

def reader_fasta(path,format_out=1):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()
	Data={}
	for index,k in enumerate(data):
		if index%2!=0:
			if format_out==1:
				Data[data[index-1].strip()[1:]]=[k]
			else:
				Data[k.strip()]=data[index-1].strip()
	return Data

def reader(path):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()
	Data=[]
	for i in data:
		if "tags" not in i:
			Data+=[i.strip().split('\t')]

	Dic={}
	for k in Data:
		if int(k[-1])>=3:
			Dic[k[0]]=float(k[-2])
	return Dic

# Read Tables
DataL1=reader("MPRA_TFBSs_1-byInsert.tsv")
DataL2=reader("MPRA_TFBSs_1-byInsert.tsv")
DataL3=reader("MPRA_TFBSs_1-byInsert.tsv")

# Extract recovered sequences
seqs=set(list(DataL1.keys())+list(DataL2.keys())+list(DataL3.keys()))
scores1=[];scores2=[];scores3=[];
for k in seqs:
	if k in DataL1.keys() and k in DataL2.keys() and k in DataL3.keys():
		scores1+=[math.log(DataL1[k],2)]
		scores2+=[math.log(DataL2[k],2)]
		scores3+=[math.log(DataL3[k],2)]		

# Correlations
print(pearsonr(scores1,scores2))
print(pearsonr(scores1,scores3))
print(pearsonr(scores2,scores3))

SeqsNames=reader_fasta("../../MPRA_design_TF_asymmetries.txt",0)

# Separate by Construct
ScoresD={}
SeqsL=reader_fasta("../new_strand.fa",1)
for name in SeqsL.keys():
	seq = SeqsL[name][0].strip()
	scoresL=[]
	if name in DataL1.keys():
		score1 = math.log(DataL1[name],2)
		scoresL+=[score1]
	if name in DataL2.keys():
		score2 = math.log(DataL2[name],2)
		scoresL+=[score2]
	if name in DataL3.keys():
		score3 = math.log(DataL3[name],2)
		scoresL+=[score3]
	if len(scoresL)>0:
		name_header=SeqsNames[seq]
		if ">Construct1" in name_header:
			name_to_use=name_header.split(">Construct1")[1]
		elif ">Construct2" in name_header:
			name_to_use=name_header.split(">Construct2")[1]

		if np.mean(scoresL)!=0:
			try:
				ScoresD[name_to_use]+=[np.mean(scoresL)]
			except:
				ScoresD[name_to_use]=[np.mean(scoresL)]

ax = plt.subplot(111)
scores1L=[];scores2L=[]
for i in ScoresD.keys():
	if len(ScoresD[i])==2:
		scores1L+=[ScoresD[i][0]]
		scores2L+=[ScoresD[i][1]]

print(spearmanr(scores1L,scores2L))
plt.scatter(scores1L,scores2L,s=7)
plt.xlabel("log"+"$_{2}$"+"(RNA/DNA) in Background 1")
plt.ylabel("log"+"$_{2}$"+"(RNA/DNA) in Background 2")
#plt.ylabel("log(RNA/DNA) in Construct 2")
plt.grid()
plt.legend(frameon=False)
plt.tight_layout()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.savefig("scatter_constructs_1_2.png")
plt.close()

