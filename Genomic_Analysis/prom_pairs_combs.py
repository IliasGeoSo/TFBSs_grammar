import re,os,sys,glob,math,pdb
from scipy.stats import binom_test
from adjustText import adjust_text

def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[] 
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

files=glob.glob("beds_promoters/*bed.plus")
#files2=files[:7]
files2=[files[int(sys.argv[1])-1]]
datafile=open("outs_pairs_promoters2/motif_pairs_"+str(int(sys.argv[1])-1)+".txt","w")
for one in files2:
        print(one)
        DataL=reader(one)+reader(one.split("plus")[0]+"minus")
        TF1=DataL[0][9]
        Dic={};Dic2={};
        for k in DataL:
                prom_id=k[10]+":"+k[11]+"-"+k[12]+"_"+k[-1]
                pos1=int(k[1])-int(k[11])
                if k[-1]=="+":
                        dist=int(k[12])-int(k[1])
                elif k[-1]=="-":
                        dist=int(k[1])-int(k[11])
                try:
                        Dic[prom_id]+=[dist]
                except:
                        Dic[prom_id]=[dist]

    
        for two in files:
                DataL2=reader(two)+reader(two.split("plus")[0]+"minus")
                TF2=DataL2[0][9]
                for k2 in DataL2:
                        prom_id=k2[10]+":"+k2[11]+"-"+k2[12]+"_"+k2[-1]
                        pos1=int(k2[1])-int(k2[11])
                        if k2[-1]=="+":
                                dist=int(k2[12])-int(k2[1])
                        elif k2[-1]=="-":
                                dist=int(k2[1])-int(k2[11])
                        try:
                                Dic2[prom_id]+=[dist]
                        except:
                                Dic2[prom_id]=[dist]

                promotersL=set(Dic.keys()+Dic2.keys())
                TF1_first=0;TF2_first=0;
                for prom in promotersL:
                        try:
                                first_set=Dic[prom]
                                second_set=Dic2[prom]
                                for inst1 in first_set:
                                        for inst2 in second_set:
                                                if inst1>inst2:
                                                        TF1_first+=1
                                                elif inst1<inst2:
                                                        TF2_first+=1
                        except:
                                pass
                sign=min(binom_test(TF1_first,TF1_first+TF2_first,0.5)*len(files)*len(files),1)

                datafile.write(TF1+'\t'+TF2+'\t'+str(TF1_first)+'\t'+str(TF2_first)+'\t'+str(sign)+'\n')
datafile.close()
