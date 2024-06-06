

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

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
plt.rcParams['axes.labelsize'] = 12
dpi_set = 72

#: G7-7501, CJ-5681, B8-5546, + CCLE/8b0b data (CJ-4923)

df_2ZA9JO = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/ab8e41e9-eff9-404d-959d-c0a6651c8480.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_BQ7050 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/16304454-a871-434f-a7c3-16d3d30b9c6a.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_CJ5681 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/2faab562-906d-4b40-bff8-2fca230c10b0.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_G77501 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/d326ba7a-991d-4d4e-ae96-1277b94fcf46.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_J78537 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/b75d54ad-33b4-45a5-956a-9a2a49c8280c.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_6DAA2E = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/6a1754ca-b4fb-408a-b2af-052ea1cd895b.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_B05705 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/afac506e-e65f-4342-9b9c-12bac2917e43.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_B85546 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/7acbc8b8-149f-4709-8be3-b47cfb1f4f78.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_BP4756 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/fc6bb84b-7b5c-4ed2-9643-8bc5bbd1402c.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_AK3444 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/c118cd4c-ead1-4bf2-8a15-dedb89fc1b98.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_BP4983 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/167e1c3d-8d85-4d05-a8f5-2f6c98da4206.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_B04818 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/2d20b2fc-fb90-49f9-99d6-1733e6da29b2.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_B05110 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/f37b1304-a096-4712-b70d-197ce3d5887f.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_B85552 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/c0aa6bcd-3dc1-476c-a2e0-7a66bcc0703c.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_CZ4863 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/31e38e13-b825-4445-99a1-d69ed4e52adf.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_G6A8L7 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/464d68d7-b51b-43a2-8a95-4b8947f28c9d.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_G6A8L8 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/92ce0913-85e6-48a2-9c55-e51489730b16.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_T7A92I = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/d2ed8261-28ac-4957-b5ea-fecdb0b73cee.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_BP4801 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/b7d93ebf-c34a-42a8-aeec-05b78f2b5495.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_CJ4887 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/dcd98863-4e58-4457-b083-dfda3c51ee33.methylation_array.sesame.level3betas.txt", sep="\t", header=None)
df_BP5170 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/methylation_subset/d73ba7a5-aa0f-49bc-b512-e0f6fc58335a.methylation_array.sesame.level3betas.txt", sep="\t", header=None)

df_ref = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/CCLE/tRCC4_high_TCGA-J7-8537-01_low_TCGA-J7-8537-01_samples.txt", sep="\t")
cpg_promoter_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/full_TCGA/chrX_non_escaping_CpG_promoter_islands.csv")
islands_of_interest = cpg_promoter_df['Composite Element REF'].tolist()

df_ref = df_ref[df_ref['chromosome']=="chrX"]

#df_ref = df_ref[df_ref['feature_type']=="Island"]
#df_ref = df_ref[df_ref['composite_element_ref'].isin(islands_of_interest)]
comp_ref = df_ref['composite_element_ref'].tolist()
pos_ref = df_ref['start'].tolist()

def generate_rows(df_temp):
    ref_temp = df_temp[0].tolist()
    val_temp = df_temp[1].tolist()
    reference_dict = dict(zip(ref_temp,val_temp))
    new_vals = []
    for a in comp_ref:
        new_vals.append(reference_dict[a])
    return new_vals

def generate_plot(value, name):
    
    try:
        df1 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/TCGA/unpaired_RNAseq_plots/TRUNCATED_withoutVline_chrX_RNA_unpaired_on_germline_hets_unphased_TCGA-" + name + ".log_AF_vs_position.csv")
    except:
        try:
            df1 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/TCGA/unpaired_RNAseq_ccRCCs_plots/TRUNCATED_withoutVline_chrX_Tophat2_RNA_unpaired_on_germline_hets_unphased_uq_mapped_TCGA-" + name + ".log_AF_vs_position.csv")
        except:
            df1 = False
    
    #Generate .csvs for ccRCC

    try:
        df1 = df1[df1['position']<49028726]
        df2 = df1[df1['MAF']<0.2]
        df1 = df1[df1['MAF']>=0.2]
        df1['Sample'] = len(df1.index.tolist())*[name]
        pos_list = df1['gene_start'].tolist()
        second_pos_list = df2['gene_start'].tolist()
    except:
        pos_list = []
        second_pos_list = []

    def generate_rows(df_temp):
        ref_temp = df_temp[0].tolist()
        val_temp = df_temp[1].tolist()
        reference_dict = dict(zip(ref_temp,val_temp))
        new_vals = []
        for a in comp_ref:
            new_vals.append(reference_dict[a])
        return new_vals

    xs = pos_ref
    ys = generate_rows(value)
    cvals = []

    #Check promoters, not just random CpGs

    for q in xs:
        x=0
        for w in pos_list:
            if q>(w-10000) and q<(w):
                x=1
        if x == 1:
            cvals.append("red")
        else:
            y=0
            for w in second_pos_list:
                if q>(w-10000) and q<(w):
                    y=1
            if y == 1:
                cvals.append("cyan")   
            else:
                cvals.append("black")

    fig = plt.figure(figsize=(5,5))
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    ax.scatter(xs,ys,c=cvals, alpha=1, s=4)

    ax.set_facecolor("white")
    ax.grid(False)
    ax.tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4)
    ax.tick_params(axis='y', which='major', left=True, labelsize=10, size=4)

    ax.spines['bottom'].set_color('0')
    ax.spines['left'].set_color('0')
    ax.spines['top'].set_color('1')
    ax.spines['right'].set_color('1')

    ax.set_xlabel("hg38 chrX Position", labelpad=10)
    ax.set_ylabel("Methylation at Probe (β)")

    for label in ax.get_xticklabels():
        label.set_ha("center")
        label.set_rotation(0)

    fig.tight_layout()

    dpi_set = 72
    plt.tick_params(bottom='on', left='on')

    ax.set_ylim([-0.01, 1.01])

    len_value = 156040895
    ax.set_xlim([-500000, len_value+500000])
    ax.axvline(49028726, lw=0.5, c="black", ls="--")

    fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/methylation_new/" + name + "_methylation_plot.pdf", dpi=dpi_set, bbox_inches = 'tight')    

    return df1

dfa = generate_plot(df_CJ5681, "CJ-5681")
dfb = generate_plot(df_2ZA9JO, "2Z-A9JO")
dfc = generate_plot(df_BQ7050, "BQ-7050")
dfd = generate_plot(df_G77501, "G7-7501")
dfe = generate_plot(df_J78537, "J7-8537")

out_df = pd.concat([dfa,dfb,dfc,dfd,dfe])
out_df.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/methylation_new/MAFgeq02_genes_in_reactivated_tRCCs.csv", index=False)

"""
generate_plot(df_6DAA2E, "6D-AA2E")
generate_plot(df_B05705, "B0-5705")
generate_plot(df_B85546, "B8-5546")
generate_plot(df_BP4756, "BP-4756")

generate_plot(df_AK3444, "AK-3444")
generate_plot(df_BP4983, "BP-4983")
generate_plot(df_B04818, "B0-4818")
generate_plot(df_B05110, "B0-5110")
generate_plot(df_B85552, "B8-5552")
generate_plot(df_CZ4863, "CZ-4863")
generate_plot(df_G6A8L7, "G6-A8L7")
generate_plot(df_G6A8L8, "G6-A8L8")
generate_plot(df_T7A92I, "T7-A92I")

generate_plot(df_BP4801, "BP-4801")
generate_plot(df_CJ4887, "CJ-4887")
generate_plot(df_BP5170, "BP-5170")
"""
