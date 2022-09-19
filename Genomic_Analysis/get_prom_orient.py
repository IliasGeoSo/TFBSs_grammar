import re,os,sys,glob

files=glob.glob("beds/*")
for one in files:
        print(one)
        os.system("bedtools intersect -a " + one + " -b gencode/gencode.v39.bed.protein_coding.promoters.bed.plus -wb > beds_promoters/"+one.split("/")[-1]+".plus")
        os.system("bedtools intersect -a " + one + " -b gencode/gencode.v39.bed.protein_coding.promoters.bed.minus -wb > beds_promoters/"+one.split("/")[-1]+".minus") 
~                                                                                                                                                                         
