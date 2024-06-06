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
sns.set(rc={'figure.figsize':(7,4)})
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
    df['Major_Allele'] = np.log2((2**df['LogRatio'])*np.maximum(df['AllelicRatio'], (1-df['AllelicRatio'])))
    df['Minor_Allele'] = np.log2((2**df['LogRatio'])*np.minimum(df['AllelicRatio'], (1-df['AllelicRatio'])))
    return df

df1 = generate_corrected_vals(df1)
df2 = generate_corrected_vals(df2)

df1['Adj_Start'] = df1['Position']
df2['Adj_Start'] = df2['Position']

xs1 = df1['Adj_Start'].tolist()
ys1 = df1['Major_Allele'].tolist()

xs2 = df1['Adj_Start'].tolist()
ys2 = df1['Minor_Allele'].tolist()

xs3 = df2['Adj_Start'].tolist()
ys3 = df2['Major_Allele'].tolist()

xs4 = df2['Adj_Start'].tolist()
ys4 = df2['Minor_Allele'].tolist()

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

fig, (ax1, ax2) = plt.subplots(2,1, gridspec_kw={'height_ratios': [1, 1]})
ax1.scatter(xs1,ys1, c=clist1, alpha=1, s=2)
ax1.scatter(xs2,ys2, c=clist2, alpha=1, s=2)
ax1.set_facecolor("white")
ax1.grid(False)
ax1.tick_params(axis='x', which='major', bottom=True, labelsize=15, size=4)
ax1.tick_params(axis='y', which='major', left=True, labelsize=15, size=4)
ax1.spines['bottom'].set_color('0')
ax1.spines['left'].set_color('0')

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

ax2.scatter(xs3,ys3, c=clist3, alpha=1, s=2)
ax2.scatter(xs4,ys4, c=clist4, alpha=1, s=2)
ax2.set_facecolor("white")
ax2.grid(False)
ax2.tick_params(axis='x', which='major', bottom=True, labelsize=15, size=4)
ax2.tick_params(axis='y', which='major', left=True, labelsize=15, size=4)
ax2.spines['bottom'].set_color('0')
ax2.spines['left'].set_color('0')

#patch1 = mpatches.Patch(color='black', label='Major Allele')
#patch2 = mpatches.Patch(color='orange', label='Minor Allele')

#leg = ax2.legend(handles=[patch1,patch2], bbox_to_anchor=(0.7, -0.4), fontsize=15, ncol=2, facecolor="white")
#frame = leg.get_frame()
#frame.set_edgecolor("black")
#frame.set_linewidth(0.5)

"""
ax1.set_ylim([-0.5, 4.5])
ax2.set_ylim([-0.5, 4.5])

ax1.set_yticks([0, 2, 4])
ax1.set_yticklabels([0, 2, 4])
ax2.set_yticks([0, 2, 4])
ax2.set_yticklabels([0, 2, 4])
"""

ax1.set_ylabel("Allelic logR", fontsize=15)
ax2.set_ylabel("Allelic logR", fontsize=15)
ax1.set_xlabel("hg38 chr1 Coordinate (bp)", fontsize=15)
ax2.set_xlabel("hg38 chrX Coordinate (bp)", fontsize=15)

ax2.axvline(49028726, ls="--", lw=2, c="black")
ax1.axvline(35176378, ls="--", lw=2, c="black")

fig.tight_layout()
dpi_set = 1200
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/SCNA_plot_TFE3-53.pdf", dpi=dpi_set)
