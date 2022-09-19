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
DicMotifs={}
for one in files:
	filed=glob.glob(one+"/*tsv")[0]
	DataL=reader_meme(filed)
	try:
		if "hoco_ahr" in one:
			TF=DataL[0][0].upper()
		else:
			TF=DataL[0][1].upper()

		TFsL+=[TF]
		for k in DataL:
			if k[2] not in DicMotifs.keys():
				DicMotifs[k[2]]={}
			if len(k)>1:
				try:	
					DicMotifs[k[2]][TF]+=[k[5]]
				except:
					DicMotifs[k[2]][TF]=[k[5]]
				#	DicMotifs[TF][k[2].replace(":","").replace("]","").replace("[","")]+=[k[5]]
				#except:
				#	DicMotifs[TF][k[2].replace(":","").replace("]","").replace("[","")]=[k[5]]
	except:
		pass
		
exp_with_without=[]
for TF in TFsL:	
	#with_=[]
	without_=[]
'''
	for line in Data.keys():
		found=DicMotifs[line][TF]
		if len(found)>0:
			with_+=[Data[line]]
		if len(found)==0:
			without_+=[Data[line]]
	if ttest_ind(without_,with_)[-1]<0.05:
		exp_with_without+=[TF]
'''
TFsL=list(set(TFsL))
Data=reader("ENCFF093PXS_hepg2.tsv")
seqs=reader_fasta("ENCFF245LAC_hepg2.fasta")
with_plusL=[];with_minusL=[];
with_plusL_two=[];with_minusL_two=[];
counts=0
datafile=open("Vikram_library_T_NT.txt","w")
datafile.write("motif"+'\t'+'name'+'\t'+'p_value'+'\t'+'Delta'+'\n')
TFs_to_use=["CTCF","NR2F2","SP1","HNF4A","YY1","FOSL1::JUN","FOXA1","TFAP2C","CREB1","PPARA::RXRA","RXRA","GABPA","HNF1A","ONECUT1","AHR::ARNT","XBP1"]
namesL=["CTCF","NR2F2","SP1","HNF4A","YY1","FOSL1::JUN","FOXA1","TFAP2C","CREB1","PPARA::RXRA","RXRA","GABPA","HNF1A","ONECUT1","AHR","XBP1","REST"]
TFs_to_use=["CTCF","NR2F2","SP1","HNF4A","YY1","FOSL1::JUN","FOXA1","TFAP2C","CREB1","PPARA::RXRA","RXRA","GABPA","HNF1A","ONECUT1","AHR_HUMAN.H11MO.0.B","XBP1","REST"]

for TF in TFsL:
	if TF.upper() in TFs_to_use:
		count=0
		with_plus=[];with_minus=[];
		with_plus_two=[];with_minus_two=[];
		datafile_=open("pics/"+TF+"_T_NT.txt","w")
		for line in Data.keys():
			#if "Reversed" not in line[2]:
				try:
					found=DicMotifs[line][TF]	
					plus=found.count("+")
					minus=found.count("-")
					if plus>0 and minus==0:
						count+=1
						with_plus+=[float(Data[line])]
						if plus>1: #and minus==0:
							with_plus_two+=[float(Data[line])]
					if minus>0 and plus==0:
						count+=1
						with_minus+=[float(Data[line])]
						if minus>1: #and plus==0:
							 with_minus_two+=[float(Data[line])]
				except:	
					pass
		if count>100:
			counts+=1
			with_plusL+=[with_plus]
			with_minusL+=[with_minus]
			with_plusL_two+=[with_plus_two]
			with_minusL_two+=[with_minus_two]
			#if ttest_ind(with_plus_two,with_minus_two)[-1]<0.05:
			#	print(TF,"here",len(with_plus_two),len(with_minus_two),ttest_ind(with_plus_two,with_minus_two))

			if ttest_ind(with_plus,with_minus)[-1]<0.05:
				print(TF,"here",len(with_plus),len(with_minus),ttest_ind(with_plus,with_minus))
		
				for vv in with_plus:
					datafile_.write(str(vv)+'\t'+"Non-Template"+'\n')
				for vv2 in with_minus:
					datafile_.write(str(vv2)+'\t'+"Template"+'\n')
				plt.boxplot([with_plus,with_minus])
				plt.savefig("pics/"+TF+"_TNT.png")
				plt.close()
			datafile_.close()

			#datafile.write(TF+'\t'+TF+'\t'+str(ttest_ind(with_plus,with_minus)[-1])+'\t'+str(abs(np.mean(with_plus)-np.mean(with_minus)))+'\n')
			datafile.write(TF+'\t'+TF+'\t'+str(ttest_ind(with_plus,with_minus)[-1])+'\t'+str(abs(np.mean(with_plus)-np.mean(with_minus))/min(math.log(abs(np.mean(with_plus)),abs(np.mean(with_minus))),10))+'\n')
datafile.close()

print(counts)

diffL=[];p_valL=[];
for k in range(len(with_plusL)):
	if with_plusL+with_minusL!=[]:
		#if mannwhitneyu(with_plusL[k],with_minusL[k])[-1]*len(exp_with_without)<0.005:
			try:
				#print(TFsL[k],with_plusL[k],with_minusL[k])
				print(np.mean(with_plusL[k]),np.mean(with_minusL[k],ttest_ind(with_plusL[k],with_minusL[k])))
				p_valL+=[100]
			except:
				pass
		#else:
		#	p_valL+=[-math.log(min(1,mannwhitneyu(with_plusL[k],with_minusL[k])[-1]*len(TFsL)))]
		#diffL+=[np.mean(with_plusL[k])-np.mean(with_minusL[k])]

#plt.scatter(diffL,p_valL,s=4,color="blue",alpha = 0.5)

'''
diffL_two=[];p_valL_two=[];
for k in range(len(with_plusL)):
        if mannwhitneyu(with_plusL_two[k],with_minusL_two[k])[-1]*len(TFsL)<10**-100:
                print(TFsL[k])
                p_valL_two+=[100]
        else:
                p_valL_two+=[-math.log(min(1,mannwhitneyu(with_plusL_two[k],with_minusL_two[k])[-1]*len(TFsL)))]
        diffL_two+=[np.mean(with_plusL_two[k])-np.mean(with_minusL_two[k])]

plt.scatter(diffL_two,p_valL_two,s=4,color="blue",alpha = 0.5)
'''
