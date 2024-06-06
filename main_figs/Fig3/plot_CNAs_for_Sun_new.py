import pandas as pd
import collections
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.ticker as mtick
from matplotlib import ticker
import seaborn as sns
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from statannot import add_stat_annotation
from scipy.stats import fisher_exact

plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(7,8)})
sns.set(font_scale=1)
plt.rcParams.update({'font.size': 15})

df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/CNIs/53_Tumor_cluster2.titan.txt", sep="\t")
df1 = df[df['Chr'].isin([str(1),1])]
df2 = df[df['Chr'].isin(["X"])]

"""
#Corrected_logR
bp_per_row = 50000

def generate_corrected_vals(df):
    #df = df[['Start','logR_Copy_Number']]
    df = df.groupby('Start').mean().reset_index()
    min_val = min(df['Start'].tolist())
    df['Adj_Start'] = df['Start'] - min_val
    df['group'] = (df['Adj_Start'] // bp_per_row) + 1
    df = df.groupby('group').apply(lambda x: x.mean())
    return df
    #df_len = len(df.index.tolist())
    #df['New_Start'] = np.arange(bp_per_row/2, df_len*bp_per_row+bp_per_row/2, bp_per_row)
    #df.to_csv("/Users/ananthansadagopan/Downloads/out.csv", index=False)
"""

def generate_corrected_vals(df):
    #df['Major_Allele'] = df['logR_Copy_Number']*np.maximum(df['Corrected_Ratio'], (1-df['Corrected_Ratio']))
    #df['Minor_Allele'] = df['logR_Copy_Number']*np.minimum(df['Corrected_Ratio'], (1-df['Corrected_Ratio']))
    df['Major_Allele'] = df['LogRatio']*np.maximum(df['AllelicRatio'], (1-df['AllelicRatio']))
    df['Minor_Allele'] = df['LogRatio']*np.minimum(df['AllelicRatio'], (1-df['AllelicRatio']))
    return df

df1 = generate_corrected_vals(df1)
df2 = generate_corrected_vals(df2)

df1['Adj_Start'] = df1['Position']
df2['Adj_Start'] = df2['Position']

xs1 = df1['Adj_Start'].tolist()
ys1 = np.maximum(df1['AllelicRatio'], (1-df1['AllelicRatio'])).tolist()

xs2 = df1['Adj_Start'].tolist()
ys2 = np.minimum(df1['AllelicRatio'], (1-df1['AllelicRatio'])).tolist()

xs3 = df2['Adj_Start'].tolist()
ys3 = np.maximum(df2['AllelicRatio'], (1-df2['AllelicRatio'])).tolist()

xs4 = df2['Adj_Start'].tolist()
ys4 = np.minimum(df2['AllelicRatio'], (1-df2['AllelicRatio'])).tolist()

ys5 = df1['LogRatio'].tolist()
ys6 = df2['LogRatio'].tolist()

clist1 = []
q=0
while q<len(xs1):
    if xs1[q]>35176378:
        clist1.append("gainsboro")
    else:
        clist1.append("black")
    q=q+1

clist2 = []
q=0
while q<len(xs2):
    if xs2[q]>35176378:
        clist2.append("gainsboro")
    else:
        clist2.append("orange")
    q=q+1

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, gridspec_kw={'height_ratios': [1, 1, 1, 1]})
ax1.scatter(xs1,ys1, c=clist1, alpha=1, s=2)
ax1.scatter(xs2,ys2, c=clist2, alpha=1, s=2)
ax2.scatter(xs1,ys5, c="gainsboro", alpha=1, s=2)

clist3 = []
q=0
while q<len(xs3):
    if xs3[q]<49028726:
        clist3.append("gainsboro")
    else:
        clist3.append("black")
    q=q+1

clist4 = []
q=0
while q<len(xs4):
    if xs4[q]<49028726:
        clist4.append("gainsboro")
    else:
        clist4.append("orange")
    q=q+1

ax3.scatter(xs3,ys3, c=clist3, alpha=1, s=2)
ax3.scatter(xs4,ys4, c=clist4, alpha=1, s=2)
ax4.scatter(xs3,ys6, c="gainsboro", alpha=1, s=2)

ax1.set_facecolor("white")
ax1.grid(False)
ax1.tick_params(axis='x', which='major', bottom=True, labelsize=15, size=4)
ax1.tick_params(axis='y', which='major', left=True, labelsize=15, size=4)
ax1.spines['bottom'].set_color('0')
ax1.spines['left'].set_color('0')

ax2.set_facecolor("white")
ax2.grid(False)
ax2.tick_params(axis='x', which='major', bottom=True, labelsize=15, size=4)
ax2.tick_params(axis='y', which='major', left=True, labelsize=15, size=4)
ax2.spines['bottom'].set_color('0')
ax2.spines['left'].set_color('0')

ax3.set_facecolor("white")
ax3.grid(False)
ax3.tick_params(axis='x', which='major', bottom=True, labelsize=15, size=4)
ax3.tick_params(axis='y', which='major', left=True, labelsize=15, size=4)
ax3.spines['bottom'].set_color('0')
ax3.spines['left'].set_color('0')

ax4.set_facecolor("white")
ax4.grid(False)
ax4.tick_params(axis='x', which='major', bottom=True, labelsize=15, size=4)
ax4.tick_params(axis='y', which='major', left=True, labelsize=15, size=4)
ax4.spines['bottom'].set_color('0')
ax4.spines['left'].set_color('0')

ax1.set_ylabel("Allelic Ratio", fontsize=15)
ax3.set_ylabel("Allelic Ratio", fontsize=15)
ax2.set_ylabel("log$_2$(T/N Ratio)", fontsize=15)
ax4.set_ylabel("log$_2$(T/N Ratio)", fontsize=15)
ax1.set_xlabel("hg38 chr1 Coordinate (bp)", fontsize=15)
ax3.set_xlabel("hg38 chrX Coordinate (bp)", fontsize=15)
ax2.set_xlabel("hg38 chr1 Coordinate (bp)", fontsize=15)
ax4.set_xlabel("hg38 chrX Coordinate (bp)", fontsize=15)

ax3.axvline(49028726, ls="--", lw=2, c="black")
ax1.axvline(35176378, ls="--", lw=2, c="black")
ax4.axvline(49028726, ls="--", lw=2, c="black")
ax2.axvline(35176378, ls="--", lw=2, c="black")

ax1.set_ylim([-0.1, 1.1])
ax3.set_ylim([-0.1, 1.1])
ax1.set_yticks([0, 0.5, 1])
ax1.set_yticklabels([0, 0.5, 1])
ax3.set_yticks([0, 0.5, 1])
ax3.set_yticklabels([0, 0.5, 1])

ax2.set_ylim([-1.1, 1.1])
ax4.set_ylim([-1.1, 1.1])
ax2.set_yticks([-1, 0, 1])
ax2.set_yticklabels([-1, 0, 1])
ax4.set_yticks([-1, 0, 1])
ax4.set_yticklabels([-1, 0, 1])

fig.tight_layout()
dpi_set = 1200
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/SCNA_plot_TFE3-53_logratios_split.pdf", dpi=dpi_set)
