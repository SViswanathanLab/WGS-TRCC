import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors
import matplotlib.patches as mpatches
import collections
import matplotlib.gridspec as gridspec
from scipy.stats import zscore
from itertools import chain
from statannot import add_stat_annotation

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
plt.rcParams['axes.labelsize'] = 12
dpi_set = 72

sns.set(rc={'figure.figsize':(3,5)})
sns.set(font_scale=1)
plt.rcParams.update({'font.size': 12})

nice_palette = ['#9D2727', '#CFB997', '#E79898', '#C445A1', '#DDDD00', '#36E2EE', '#008BB5', '#DACAFF', '#FFCAFF']

df_RNA = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/TCGA/TCGA_all_tRCCs_RNAseq.csv")
all_RNA_cols = df_RNA.columns.tolist()

new_RNA_cols = []
for a in all_RNA_cols:
    new_RNA_cols.append(a.split("|")[0])

df_RNA.columns = new_RNA_cols

reactivated = ['TCGA-2Z-A9JO-01A-11R-A42S-07','TCGA-BQ-7050-01A-11R-1965-07','TCGA-G7-7501-01A-11R-2204-07','TCGA-J7-8537-01A-11R-2404-07','TCGA-CJ-5681-01A-11R-1541-07']
df1_RNA = df_RNA[df_RNA['sample'].isin(reactivated)]
df2_RNA = df_RNA[~(df_RNA['sample'].isin(reactivated))]

all_reactivated_samp_names = df1_RNA['sample'].tolist()

dfvals = ['2Z-A9JO','BQ-7050','CJ-5681','G7-7501','J7-8537']

data_df = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/methylation_new/MAFgeq02_genes_in_reactivated_tRCCs_all.xlsx")
data_df_sub = data_df.drop_duplicates(subset='gene_start',keep="first")
uq_starts = [int(x) for x in data_df_sub['gene_start'].tolist()]
uq_genes = list(set(data_df_sub['gene'].tolist()))
uq_genes.remove('ARSL')
uq_genes.append('ARSE')

norm_val_R = []
norm_val_NR = []

q=0
while q<len(uq_genes):

    data_tempdf = data_df[data_df['gene']==uq_genes[q]]
    reactivated_samps = data_tempdf['Sample'].tolist()

    valid_samp_names = []
    for a in all_reactivated_samp_names:
        temp_val = split_advanced(split_advanced(a, "-", 3)[0], "-", 1)[1]
        if temp_val in reactivated_samps:
            valid_samp_names.append(a)
    
    dfR_RNA = df1_RNA[df1_RNA['sample'].isin(valid_samp_names)]

    try:
        dfR_RNA_exp = dfR_RNA[uq_genes[q]].tolist()
        dfNR_RNA_exp = df2_RNA[uq_genes[q]].tolist()
    except:
        q=q+1
        continue

    mean1 = np.mean(dfNR_RNA_exp)

    for a in dfR_RNA_exp:
        if a/mean1 >= 2.6:
            print(a/mean1)
            print(uq_genes[q])
        norm_val_R.append(a/mean1)
    
    for a in dfNR_RNA_exp:
        norm_val_NR.append(a/mean1)

    q=q+1

master_df = pd.DataFrame([norm_val_R, norm_val_NR]).T
column_list = ['Reactivated', 'Non-Reactivated']
master_df.columns = column_list
print(master_df)

fig = plt.figure()
gs = gridspec.GridSpec(1,1)
ax_stripplot = fig.add_subplot(gs[0])

nice_palette = ['#9D2727', '#CFB997', '#E79898', '#C445A1', '#DDDD00', '#36E2EE', '#008BB5', '#DACAFF', '#FFCAFF']
sns.boxplot(data=master_df, ax=ax_stripplot, fliersize=10, showfliers=False, palette=nice_palette)
sns.stripplot(data=master_df, color="black", s=4, ax=ax_stripplot, jitter=0.2, edgecolor="black", alpha=0.9, zorder=5, linewidth=0.15)

ax_stripplot.set(ylabel="Normalized Gene Expression")
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

add_stat_annotation(ax_stripplot, data=master_df,
                    box_pairs=[(column_list[0], column_list[1])],
                    test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=1/20, line_offset=0.01)


ax_stripplot.grid(False)
ax_stripplot.set_facecolor("white")

ax_stripplot.spines['bottom'].set_color('0')
ax_stripplot.spines['left'].set_color('0')

plt.tick_params(bottom='on', left='on')

fig.tight_layout()

dpi_set = 1200

fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/TCGA/TCGA_reactivated_vs_nonreactivated_normalized_gene_expression.pdf", dpi=dpi_set)
    

    

