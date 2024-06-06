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

plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(7,10)})
sns.set(font_scale=1)
plt.rcParams.update({'font.size': 13})

weird_fusions = ['TCGA-CJ-5681_Tumor','tRCC_72_T', 'tRCC_03_T', 'tRCC_23_T', 'tRCC_32_T', 'tRCC_37_T', 'tRCC_51_T', '04_Tumor', '18_Tumor', 'TCGA-BQ-5887_Tumor']
chr_name = [19, 17, 22, 17, 1, 17, 17, 17, 1, 1]
pos = [6413102, 81976807, 20495913, 81976807, 156750610, 81976807, 81976807, 50184101, 156750610, 156750610]

values = os.listdir("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/CNIs/")
values_to_iterate = [x for x in values if "titan.ichor.cna.txt" in x and ".pdf" not in x]

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

    df_seg = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/CNIs/" + temp_value[:-8] + ".seg.txt", sep="\t")
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

    xs1 = df1['Adj_Start'].tolist()
    ys1 = np.maximum(df1['AllelicRatio'], (1-df1['AllelicRatio'])).tolist()

    xs2 = df1['Adj_Start'].tolist()
    ys2 = np.minimum(df1['AllelicRatio'], (1-df1['AllelicRatio'])).tolist()

    xs3 = df2['Adj_Start'].tolist()

    ys5 = 2**np.array(df1['LogRatio'].tolist())
    ys6 = 2**np.array(df2['LogRatio'].tolist())

    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5,1, gridspec_kw={'height_ratios': [1, 1, 1, 1, 1]})

    copyneutral_starts = []
    copyneutral_ends = []
    start_vals = df_seg1['Start'].tolist()
    end_vals = df_seg1['End'].tolist()
    alleleA_CN = df_seg1['Corrected_MajorCN'].tolist()
    alleleB_CN = df_seg1['Corrected_MinorCN'].tolist()

    qx=0
    while qx<len(start_vals):
        if float(alleleA_CN[qx]) == float(alleleB_CN[qx]):
            ax3.plot((float(start_vals[qx]), float(end_vals[qx])), (float(alleleA_CN[qx]), float(alleleA_CN[qx])), marker=None, linestyle="-", lw=5, c="gainsboro", alpha=1)
            ax3.plot((float(start_vals[qx]), float(end_vals[qx])), (float(alleleB_CN[qx]), float(alleleB_CN[qx])), marker=None, linestyle="-", lw=5, c="gainsboro", alpha=0.6)
            copyneutral_starts.append(float(start_vals[qx]))
            copyneutral_ends.append(float(end_vals[qx]))
        else:
            ax3.plot((float(start_vals[qx]), float(end_vals[qx])), (float(alleleA_CN[qx]), float(alleleA_CN[qx])), marker=None, linestyle="-", lw=5, c="black", alpha=1)
            ax3.plot((float(start_vals[qx]), float(end_vals[qx])), (float(alleleB_CN[qx]), float(alleleB_CN[qx])), marker=None, linestyle="-", lw=5, c="orange", alpha=0.6)
        qx=qx+1

    if "tRCC_19_T" in temp_value:
        copyneutral_starts = copyneutral_starts[:-1]
        copyneutral_starts.append(default_pos)

    valid_pos = []
    all_positions = df1['Adj_Start'].tolist()
    w=0
    while w<len(copyneutral_starts):
        for z in all_positions:
            if z >= int(copyneutral_starts[w]) and z <= int(copyneutral_ends[w]):
                valid_pos.append(z)
        w=w+1

    df1_cn = df1[df1['Adj_Start'].isin(valid_pos)]
    new_pos = df1_cn['Adj_Start'].tolist()
    ys1_cn = np.maximum(df1_cn['AllelicRatio'], (1-df1_cn['AllelicRatio'])).tolist()
    ys2_cn = np.minimum(df1_cn['AllelicRatio'], (1-df1_cn['AllelicRatio'])).tolist()

    df1_inverse = df1[~(df1['Adj_Start'].isin(valid_pos))]
    org_pos = df1_inverse['Adj_Start'].tolist()
    org_pos_cn = np.minimum(df1_inverse['AllelicRatio'], (1-df1_inverse['AllelicRatio'])).tolist()

    #ax1.scatter(xs1,ys1, c="black", alpha=1, s=2)
    ax1.scatter(org_pos,org_pos_cn, c="orange", s=2, alpha=1)
    #ax1.scatter(new_pos,ys1_cn, c="gainsboro", alpha=1, s=2.01)
    ax1.scatter(new_pos,ys2_cn, c="gainsboro", s=2, alpha=1)
    ax2.scatter(xs1,ys5, c="gainsboro", s=2, alpha=1)

    ax4.scatter(xs3,ys6, c="black", alpha=1, s=2)

    start_vals = df_seg2['Start'].tolist()
    end_vals = df_seg2['End'].tolist()
    alleleA_CN = df_seg2['Corrected_Copy_Number'].tolist()

    qx=0
    while qx<len(start_vals):
        ax5.plot((float(start_vals[qx]), float(end_vals[qx])), (float(alleleA_CN[qx]), float(alleleA_CN[qx])), marker=None, linestyle="-", lw=5, c="black", alpha=1)
        qx=qx+1

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

    ax5.set_facecolor("white")
    ax5.grid(False)
    ax5.tick_params(axis='x', which='major', bottom=True, labelsize=15, size=4)
    ax5.tick_params(axis='y', which='major', left=True, labelsize=15, size=4)
    ax5.spines['bottom'].set_color('0')
    ax5.spines['left'].set_color('0')

    ax1.set_ylim([-0.1, 0.6])
    ax1.set_yticks([0, 0.25, 0.5])
    ax1.set_yticklabels([0, 0.25, 0.5])

    ax2.set_ylim([0, 2])
    ax2.set_yticks([0, 1, 2])
    ax2.set_yticklabels([0, 1, 2])

    ax4.set_ylim([0, 2])
    ax4.set_yticks([0, 1, 2])
    ax4.set_yticklabels([0, 1, 2])

    ax3.set_ylim([-0.5, 6.1])
    ax5.set_ylim([-0.5, 6.1])
    ax3.set_yticks([0, 2, 4, 6])
    ax3.set_yticklabels([0, 2, 4, 6])
    ax5.set_yticks([0, 2, 4, 6])
    ax5.set_yticklabels([0, 2, 4, 6])

    ax1.set_ylabel("Minor\nAllele Fraction", fontsize=15)
    ax2.set_ylabel("Copy Ratio", fontsize=15)
    ax4.set_ylabel("Copy Ratio", fontsize=15)
    ax3.set_ylabel("Allelic\nCopy Number", fontsize=15)
    ax5.set_ylabel("Copy Number", fontsize=15)

    ax1.set_xlabel("hg38 chr" + str(default_chr) + " Coordinate (bp)", fontsize=15)
    ax2.set_xlabel("hg38 chr" + str(default_chr) + " Coordinate (bp)", fontsize=15)
    ax3.set_xlabel("hg38 chr" + str(default_chr) + " Coordinate (bp)", fontsize=15)
    ax4.set_xlabel("hg38 chrX Coordinate (bp)", fontsize=15)
    ax5.set_xlabel("hg38 chrX Coordinate (bp)", fontsize=15)

    ax1.axvline(default_pos, ls="--", lw=2, c="black")
    ax2.axvline(default_pos, ls="--", lw=2, c="black")
    ax3.axvline(default_pos, ls="--", lw=2, c="black")
    ax4.axvline(49028726, ls="--", lw=2, c="black")
    ax5.axvline(49028726, ls="--", lw=2, c="black")

    fig.tight_layout()
    plt.subplots_adjust(hspace=1.2)

    dpi_set = 300
    fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/CNIs/REVISED_" + temp_value + ".pdf", dpi=dpi_set)
    #fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/CNIs/" + temp_value + ".png", dpi=dpi_set)
