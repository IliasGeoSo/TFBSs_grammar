import re,os,sys,glob

os.system("fimo -oc hoco"+"/ --bfile background_model.hg38 --thresh 1e-4 --max-stored-scores 10000000 HOCOMOCOv11_core_HUMAN_mono_meme_format.meme ENCFF245LAC_hepg2.fasta")
