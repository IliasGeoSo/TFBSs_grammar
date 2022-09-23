import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os, glob, re
import numpy as np

plt.figure(figsize=(3, 5.8))
plt.subplots_adjust(hspace=0.20)
plt.rcParams["font.family"] = "Helvetica"
sns.set_style("white")
plt.grid(color='grey', linestyle='-', linewidth=0.5)
ax.set_ylim(-1,0.5)

filename="T_NT_AGTTAATGATTAACCAA_1"
df=pd.read_csv(filename, sep="\t", header=None)
df.columns=['distance','log','TF', "Strand"]
df['distance'] = df['distance'].apply(lambda x: x-15)

dist=df.iloc[:,0].to_numpy()
bins=[40, 80, 120, 160, 250]
bin_indices = np.digitize(dist, bins)
df['bin']=bin_indices

TF=list(set(df["TF"]))[0]

hue_plot_params = {
'data': df,
'x': 'bin',
'y': 'log',
"order": [0,1,2,3,4],
"hue": "Strand",
"hue_order" : ["Non-Template", "Template"],
"palette": {"Template":'#FFCE54',"Non-Template":'#8731C2'},
'dodge': True
}

# Add stripplot
sns.stripplot(ax=ax, **hue_plot_params, alpha=0.8,size=1.5, jitter=False)
# Add boxplot
sns.boxplot(ax=ax, **hue_plot_params, showfliers = False, showmeans=False, linewidth=0.5)

for patch in ax.artists:
  r, g, b, a = patch.get_facecolor()
  patch.set_facecolor((r, g, b, 0))

ax.set_ylabel("log2(RNA/DNA)", fontsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.set_title(TF, fontsize=6)
ax.get_legend().remove()

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.tick_params(bottom=True)
ax.set_xticklabels(["(-200,-160)", "(-160,-120)", "(-120,-80)", "(-80,-40)", "(-40,0)"])
ax.xaxis.set_tick_params(labelsize=6)
ax.set_xlabel("Distance (bp)", fontsize=6)

plt.show()
