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
import os
import statistics

#Samples to revisit: tRCC_15, 50_Tumor + 04_Tumor (has no X)

plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(7,4)})
sns.set(font_scale=1)
plt.rcParams.update({'font.size': 13})

weird_fusions = ['TCGA-CJ-5681_Tumor','tRCC_72_T', 'tRCC_03_T', 'tRCC_23_T', 'tRCC_32_T', 'tRCC_37_T', 'tRCC_51_T', '04_Tumor', '18_Tumor', 'TCGA-BQ-5887_Tumor']
chr_name = [19, 17, 22, 17, 1, 17, 17, 17, 1, 1]
pos = [6413102, 81976807, 20495913, 81976807, 156750610, 81976807, 81976807, 50184101, 156750610, 156750610]
#tRCC_45_T Fusion partner unknown

values = os.listdir("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/CNIs/")
values_to_iterate = [x for x in values if "titan.txt" in x and ".pdf" not in x and ".png" not in x]

for temp_value in values_to_iterate:

    print(temp_value)
    df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/CNIs/" + temp_value, sep="\t")

    default_chr = 1
    default_pos = 35176378

    sub_val = temp_value.split("_cluster")[0]
    if sub_val in weird_fusions:
        qq=0
        while qq<len(weird_fusions):
            if sub_val == weird_fusions[qq]:
                default_chr = chr_name[qq]
                default_pos = pos[qq]
            qq=qq+1

    df1 = df[df['Chr'].isin(["chr" + str(default_chr),str(default_chr),default_chr])]
    df2 = df[df['Chr'].isin(["chrX","X"])]

    df_seg = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/CNIs/" + temp_value[:-10] + ".segs.txt", sep="\t")
    df_seg1 = df_seg[df_seg['Chromosome'].isin(["chr" + str(default_chr),str(default_chr),default_chr])]
    df_seg2 = df_seg[df_seg['Chromosome'].isin(["chrX","X"])]

    def generate_corrected_vals(df):
        df['Major_Allele'] = df['logR_Copy_Number']*np.maximum(df['Corrected_Ratio'], (1-df['Corrected_Ratio']))
        df['Minor_Allele'] = df['logR_Copy_Number']*np.minimum(df['Corrected_Ratio'], (1-df['Corrected_Ratio']))
        return df

    df1 = generate_corrected_vals(df1)
    df2 = generate_corrected_vals(df2)

    df1['Adj_Start'] = df1['Position']
    df2['Adj_Start'] = df2['Position']

    n_overall=50

    xs1 = df1['Adj_Start'].tolist()
    ys1 = df1['Major_Allele'].tolist()

    n1 = int(len(xs1)/n_overall)

    #xs1 = [statistics.mean(xs1[i:i + n1]) for i in range(0, len(xs1), n1)]
    #ys1 = [round(statistics.mean(ys1[i:i + n1])) for i in range(0, len(ys1), n1)]

    xs2 = df1['Adj_Start'].tolist()
    ys2 = df1['Minor_Allele'].tolist()

    #xs2 = [statistics.mean(xs2[i:i + n1]) for i in range(0, len(xs2), n1)]
    #ys2 = [round(statistics.mean(ys2[i:i + n1])) for i in range(0, len(ys2), n1)]

    xs3 = df2['Adj_Start'].tolist()
    ys3 = df2['Major_Allele'].tolist()

    n2 = int(len(xs3)/n_overall)

    #xs3 = [statistics.mean(xs3[i:i + n2]) for i in range(0, len(xs3), n2)]
    #ys3 = [round(statistics.mean(ys3[i:i + n2])) for i in range(0, len(ys3), n2)]

    xs4 = df2['Adj_Start'].tolist()
    ys4 = df2['Minor_Allele'].tolist()

    #xs4 = [statistics.mean(xs4[i:i + n2]) for i in range(0, len(xs2), n2)]
    #ys4 = [round(statistics.mean(ys4[i:i + n2])) for i in range(0, len(ys2), n2)]

    fig, (ax1, ax2) = plt.subplots(2,1, gridspec_kw={'height_ratios': [1, 1]})

    #ax1.scatter(xs1,ys1, c="black", alpha=1, s=2)
    #ax1.scatter(xs2,ys2, c="orange", alpha=1, s=2)

    start_vals = df_seg1['Start_Position.bp.'].tolist()
    end_vals = df_seg1['End_Position.bp.'].tolist()
    alleleA_CN = df_seg1['Corrected_MajorCN'].tolist()
    alleleB_CN = df_seg1['Corrected_MinorCN'].tolist()

    qx=0
    while qx<len(start_vals):
        ax1.plot((float(start_vals[qx]), float(end_vals[qx])), (float(alleleA_CN[qx]), float(alleleA_CN[qx])), marker=None, linestyle="-", lw=5, c="black", alpha=1)
        ax1.plot((float(start_vals[qx]), float(end_vals[qx])), (float(alleleB_CN[qx]), float(alleleB_CN[qx])), marker=None, linestyle="-", lw=5, c="orange", alpha=0.6)
        qx=qx+1

    ax1.set_facecolor("white")
    ax1.grid(False)
    ax1.tick_params(axis='x', which='major', bottom=True, labelsize=15, size=4)
    ax1.tick_params(axis='y', which='major', left=True, labelsize=15, size=4)
    ax1.spines['bottom'].set_color('0')
    ax1.spines['left'].set_color('0')

    #ax2.scatter(xs3,ys3, c="black", alpha=1, s=2)
    #ax2.scatter(xs4,ys4, c="orange", alpha=1, s=2)

    start_vals = df_seg2['Start_Position.bp.'].tolist()
    end_vals = df_seg2['End_Position.bp.'].tolist()
    alleleA_CN = df_seg2['Corrected_MajorCN'].tolist()
    alleleB_CN = df_seg2['Corrected_MinorCN'].tolist()

    qx=0
    while qx<len(start_vals):
        ax2.plot((float(start_vals[qx]), float(end_vals[qx])), (float(alleleA_CN[qx]), float(alleleA_CN[qx])), marker=None, linestyle="-", lw=5, c="black", alpha=1)
        ax2.plot((float(start_vals[qx]), float(end_vals[qx])), (float(alleleB_CN[qx]), float(alleleB_CN[qx])), marker=None, linestyle="-", lw=5, c="orange", alpha=0.6)
        qx=qx+1

    ax2.set_facecolor("white")
    ax2.grid(False)
    ax2.tick_params(axis='x', which='major', bottom=True, labelsize=15, size=4)
    ax2.tick_params(axis='y', which='major', left=True, labelsize=15, size=4)
    ax2.spines['bottom'].set_color('0')
    ax2.spines['left'].set_color('0')

    ax1.set_ylim([-0.5, 6.1])
    ax2.set_ylim([-0.5, 6.1])
    ax1.set_yticks([0, 2, 4, 6])
    ax1.set_yticklabels([0, 2, 4, 6])
    ax2.set_yticks([0, 2, 4, 6])
    ax2.set_yticklabels([0, 2, 4, 6])

    ax1.set_ylabel("Allelic\nCopy Number", fontsize=15)
    ax2.set_ylabel("Allelic\nCopy Number", fontsize=15)
    ax1.set_xlabel("hg38 chr" + str(default_chr) + " Coordinate (bp)", fontsize=15)
    ax2.set_xlabel("hg38 chrX Coordinate (bp)", fontsize=15)

    ax2.axvline(49028726, ls="--", lw=2, c="black")
    ax1.axvline(default_pos, ls="--", lw=2, c="black")

    fig.tight_layout()
    plt.subplots_adjust(hspace=1.2)
    dpi_set = 300
    fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/CNIs/SEG_only_" + temp_value + ".pdf", dpi=dpi_set)
    #fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/CNIs/" + temp_value + ".png", dpi=dpi_set)
