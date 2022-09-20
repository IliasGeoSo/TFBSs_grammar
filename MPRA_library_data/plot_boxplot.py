import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
import pandas as pd
import os

TF = "TF"
filepath = "boxplot_txt/main/TF.txt"

plt.rcParams["font.family"] = "Helvetica"
plt.figure(figsize=(2, 10))
plt.subplots_adjust(hspace=0.25)
plt.grid(color='grey', linestyle='-', linewidth=0.5)

df=pd.read_csv(filepath, sep="\t", header=None)
df.columns=['log', 'dir']

sns.set_style("white")
order = ['Non-Template', 'Divergent', 'Convergent', 'Template']

ax = sns.boxplot(data=df,
				x='dir',
				y='log', 
				order=order, 
				showfliers = False, 
				width = 0.3,
				showmeans = False, 
				linewidth = 0.5
				)

sns.stripplot(data = df,
				x ='dir',
				y ='log', 
				order = order,
				alpha = 0.8,
				size = 1.5,
				jitter = 0,
				palette="husl"
				)

add_stat_annotation(ax, 
					data=df, 
					x='dir', 
					y='log', 
					order=order,
					box_pairs=[('Non-Template', 'Divergent'), ('Non-Template', 'Convergent'), ('Convergent', 'Divergent'), ('Convergent', 'Template'), ('Template', 'Divergent'),('Non-Template', 'Template')],
					test='t-test_ind'
					)
