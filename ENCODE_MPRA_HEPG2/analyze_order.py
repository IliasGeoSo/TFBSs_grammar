import re,os,sys,glob,math
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from scipy.stats import pearsonr,ttest_ind,kruskal
from scipy.stats import wilcoxon
import numpy as np
from Bio.Seq import Seq
from matplotlib.pyplot import figure

def reader_fasta(path):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()
	Data={}
	for index,k in enumerate(data):
		if index%2!=0:
			#Data[data[index-1].strip()[1:].replace(":","").replace("[","").replace("]","")]=k.strip()
			Data[data[index-1].strip()[1:]]=k.strip()
	return Data

def reader_fasta2(path):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()
	Data={}
	for index,k in enumerate(data):
		if index%2!=0:
			Data[k[0].strip()]=data[index-1].strip()
	return Data

def reader(path):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()
	Data=[]
	for i in data:
		if i[0]!="#":
			Data+=[i.strip().split('\t')]

	Dic={}
	for k in Data:
		if k[0]!="name":
			Dic[k[0].strip()]=float(k[1].strip())
			#Dic[k[0].strip()[1:].replace(":","").replace("[","").replace("]","")]=float(k[1].strip())
	return Dic

def reader_meme(path):
	datafile=open(path,"r")
	data=datafile.readlines()
	datafile.close()
	Data=[]
	for i in data[1:]:
		if i[0]!="#" and len(i)>1:
			if float(i.split("\t")[-3])<=10**-4:
				Data+=[i.strip().split('\t')]
	return Data


files=glob.glob("outs2/*")+["/Users/iliasgeoso/Downloads/vikrams_mpras/hoco_ahr/"]
TFsL=[]
DicMotifs={};
for one in files:
	filed=glob.glob(one+"/*tsv")[0]
	DataL=reader_meme(filed)
	try:
		if "hoco_ahr" in one:
			TF=DataL[0][0].upper()
		else:
			TF=DataL[0][1].upper()
		print(TF)
		TFsL+=[TF]
		for k in DataL:
			# start is k[3]
			if k[2] not in DicMotifs.keys():
				DicMotifs[k[2]]={}
			if len(k)>1:
				try:	
					DicMotifs[k[2]][TF]+=[k[3]]
				except:
					DicMotifs[k[2]][TF]=[k[3]]
				#	DicMotifs[TF][k[2].replace(":","").replace("]","").replace("[","")]+=[k[3]]
				#except:
				#	DicMotifs[TF][k[2].replace(":","").replace("]","").replace("[","")]=[k[3]]
	except:
		pass
		

TFsL=list(set(TFsL))
Data=reader("ENCFF093PXS_hepg2.tsv")
seqs=reader_fasta("ENCFF245LAC_hepg2.fasta")
counts=0
#datafile=open("Vikram_library_T_NT.txt","w")
#datafile.write("motif"+'\t'+'name'+'\t'+'p_value'+'\t'+'Delta'+'\n')
namesL=["CTCF","NR2F2","SP1","HNF4A","YY1","FOSL1::JUN","FOXA1","TFAP2C","CREB1","PPARA::RXRA","RXRA","GABPA","HNF1A","ONECUT1","AHR","XBP1","REST"]
TFs_to_use=["CTCF","NR2F2","SP1","HNF4A","YY1","FOSL1::JUN","FOXA1","TFAP2C","CREB1","PPARA::RXRA","RXRA","GABPA","HNF1A","ONECUT1","AHR_HUMAN.H11MO.0.B","XBP1","REST"]
AllL=[];AllpL=[];counts=0;
for TF1 in TFs_to_use:
	All=[];Allp=[];
	for TF2 in TFs_to_use:
		T1T2L=[];T2T1L=[];
		if TF1!=TF2:
			count=0
			#datafile_=open("pics2/"+TF1+"_"+TF2+".txt","w")
			for line in Data.keys():
				try:
					pos1L=DicMotifs[line][TF1]
					pos2L=DicMotifs[line][TF2]
					for pos1 in pos1L:
						for pos2 in pos2L:
							if int(pos1)<int(pos2):
								T1T2L+=[float(Data[line])]
							else:
								T2T1L+=[float(Data[line])]
				except:	
					pass
			print(mannwhitneyu(T1T2L,T2T1L))
			All+=[np.mean(T1T2L)-np.mean(T2T1L)]
			if ttest_ind(T1T2L,T2T1L)[-1]*len(TFs_to_use)*(len(TFs_to_use)-1)<0.001:
				Allp+=["***"];
				counts+=1;
			elif ttest_ind(T1T2L,T2T1L)[-1]*len(TFs_to_use)*(len(TFs_to_use)-1)<0.01:
				Allp+=["**"];
				counts+=1;
			elif ttest_ind(T1T2L,T2T1L)[-1]*len(TFs_to_use)*(len(TFs_to_use)-1)<0.05:
				Allp+=["*"];
				counts+=1;
			else:
				Allp+=[" "];

		else:
			All+=[0]
			Allp+=[" "]
	print(All)
	AllL+=[All]
	AllpL+=[Allp]

print(counts)
import seaborn as sns
import pandas as pd
df=pd.DataFrame(AllL,index=namesL,columns=namesL)
df=df.replace(np.nan, 0)
df2=pd.DataFrame(AllpL,index=namesL,columns=namesL)
g=sns.clustermap(df,annot=df2,center=0,cmap="PuOr",fmt='',linewidth=0.5)
plt.savefig("order_MPRA_vikram_hepg2.png")
plt.close()
