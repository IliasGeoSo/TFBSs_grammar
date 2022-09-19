import re,os,sys,glob

files=glob.glob("file_memes/*.meme")
motif = files[int(sys.argv[1])-1]
os.system("mkdir outs_e4/"+motif.split("/")[-1]
os.system("~/fimo_run/meme-5.0.5/src/fimo -oc outs_e4/"+motif.split("/")[-1]+"/ --bfile ../../All_motifs_strand/background_model.hg38 --thresh 1e-4 --max-stored-scores 10000000 "+motif +" hg38.fa")
