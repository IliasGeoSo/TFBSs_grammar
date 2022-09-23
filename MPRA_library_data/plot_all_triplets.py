file=["all_NR2F2_HNF1A_REST.txt"]
TF=file.split("all_")[1].split(".txt")[0]

sns.set(font="Helvetica")
sns.set_style("white")
plt.grid(color='grey', linestyle='-', linewidth=0.5)
ax.set_xlim(-1.5,1)

df=pd.read_csv(file, sep="\t", header=None)
df.columns=['log', 'dir']

group_means = df.groupby(['dir'])['log'].mean()
group_means.sort_values(inplace=True)
order = group_means.index.to_list()
pairs = list(itertools.combinations(order, 2))
pairs.sort(reverse=True)

hue_plot_params = {
  'data': df,
  'x': 'log',
  'y': 'dir',
  "order": order,
  "palette": "husl",
  'dodge': True, 
  'orient': 'h',
  }

sns.boxplot(ax=ax, 
            data=df, 
            x='log', 
            y='dir', 
            order=order, 
            orient='h', 
            showfliers = False, 
            palette='husl',
            width=0.3,
            showmeans=False,
            linewidth=0.5)

sns.stripplot(ax=ax, **hue_plot_params, alpha=0.8,size=1.5, jitter=0)

ax.set_xlabel("log2(RNA/DNA)", fontsize=8)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.title.set_text(TF.replace("_","-"))

plt.show()
