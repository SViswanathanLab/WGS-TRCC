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
import scipy

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

plt.rcParams["font.family"] = "Arial"
plt.rcParams['axes.labelsize'] = 12
dpi_set = 72
nice_palette = ['#9D2727', '#CFB997', '#E79898', '#C445A1', '#DDDD00', '#36E2EE', '#008BB5', '#DACAFF', '#FFCAFF']

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
coord_dict = dict(zip(comp_ref,pos_ref))

data_df = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/methylation_new/MAFgeq02_genes_in_reactivated_tRCCs_all.xlsx")
data_df_sub = data_df.drop_duplicates(subset='gene_start',keep="first")
#data_df = data_df[data_df['lit_annotation']!="E"]
uq_starts = [int(x) for x in data_df_sub['gene_start'].tolist()]
uq_genes = data_df_sub['gene'].tolist()

dflist = [df_2ZA9JO,df_BQ7050,df_CJ5681,df_G77501,df_J78537,df_6DAA2E,df_B05705,df_B85546,df_BP4756]
dfvals = ['2Z-A9JO','BQ-7050','CJ-5681','G7-7501','J7-8537','6D-AA2E','B0-5705','B8-5546','BP-4756']
Xa_samples = ['6D-AA2E','B0-5705','B8-5546','BP-4756']

spacer1 = 1000
spacer2 = 1000

biallelic_betas = []
monoallelic_betas = []

q=0
while q<len(dfvals):
    reads_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/TCGA/unpaired_RNAseq_plots/TRUNCATED_withoutVline_chrX_RNA_unpaired_on_germline_hets_unphased_TCGA-" + dfvals[q] + ".log_AF_vs_position.csv")
    temp_df = dflist[q]
    print(dfvals[q])
    
    altdf = reads_df[reads_df['MAF']>=0.2]
    invalid_genes = altdf['gene'].tolist()

    xpdf = reads_df[reads_df['position']<49028726]
    xpdf = xpdf[xpdf['MAF']>=0.2]
    biallelicXp_locations = list(set([x for x in xpdf['gene_start'].tolist() if x == x]))

    #Option #1
    #monoallelicdf = reads_df[reads_df['position']<49028726]

    #Option #2
    monoallelicdf = reads_df[reads_df['MAF']<0.2]
    monoallelicdf = monoallelicdf[~(monoallelicdf['gene'].isin(invalid_genes))]
    monoallelic_locations = list(set([x for x in monoallelicdf['gene_start'].tolist() if x == x]))
    monoallelic_locations = list(np.setdiff1d(monoallelic_locations,biallelicXp_locations))

    temp_dfprobes = temp_df[0].tolist()
    new_coords = []

    fig = plt.figure(figsize=(5,5))
    gs = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0])

    for a in temp_dfprobes:
        try:
            new_coords.append(coord_dict[a])
        except:
            new_coords.append(0)
        
    temp_df['coord'] = new_coords
    
    ys1 = []
    xs1 = []

    if dfvals[q] not in Xa_samples:
        for a in biallelicXp_locations:
            sub_df = temp_df[temp_df['coord']<a+spacer1]
            sub_df = sub_df[sub_df['coord']>a-spacer2]

            beta_vals = sub_df[1].tolist()
            newposvals = sub_df['coord'].tolist()

            for w in beta_vals:
                biallelic_betas.append(w)
                ys1.append(w)
            for w in newposvals:
                xs1.append(w) 
    
    ys2 = []
    xs2 = []
    
    for a in monoallelic_locations:
        sub_df = temp_df[temp_df['coord']<a+spacer1]
        sub_df = sub_df[sub_df['coord']>a-spacer2]

        beta_vals = sub_df[1].tolist()
        print(beta_vals)
        newposvals = sub_df['coord'].tolist()

        for w in beta_vals:
            monoallelic_betas.append(w)
            ys2.append(w)
        for w in newposvals:
            xs2.append(w)
    
    """
    ax.scatter(xs2, ys2, s=5, c='black')
    ax.scatter(xs1, ys1, s=5, c='red')
    
    u,p = scipy.stats.mannwhitneyu([x for x in ys1 if x==x],[x for x in ys2 if x==x])
    ax.text(0.02, 0.98, f'P = {p:.3f}', transform=ax.transAxes, fontsize=14, verticalalignment='top')
    
    ax.set_facecolor("white")
    ax.grid(False)
    ax.tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4)
    ax.tick_params(axis='y', which='major', left=True, labelsize=10, size=4)

    ax.spines['bottom'].set_color('0')
    ax.spines['left'].set_color('0')
    ax.spines['top'].set_color('1')
    ax.spines['right'].set_color('1')
    ax.set_ylabel("Methylation at Probe (β)")
    ax.set_xlabel("hg38 chrX Coordinate")

    for label in ax.get_xticklabels():
        label.set_ha("center")
        label.set_rotation(0)

    ax.set_ylim([-0.01, 1.01])
    #ax.axhline(0.4, ls="--", lw=0.5, c="black")

    fig.tight_layout()

    dpi_set = 72
    plt.tick_params(bottom='on', left='on')

    fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/methylation_new/methylation_sample_level_" + dfvals[q] + ".pdf", dpi=dpi_set, bbox_inches = 'tight')    
    """
    q=q+1

fig = plt.figure(figsize=(3,5))
gs = gridspec.GridSpec(1,1)
ax = fig.add_subplot(gs[0])

print(biallelic_betas)
print(monoallelic_betas)

master_df = pd.DataFrame([[float(x) for x in biallelic_betas if x==x],[float(x) for x in monoallelic_betas if x==x]]).T
col_names = ['MAF≥0.2\nchrXp','MAF<0.2']
master_df.columns = col_names

print(master_df)

"""
validgeqcutoff = [x for x in all_valid_values if x == x and x > 0.4]
invalidgeqcutoff = [x for x in all_invalid_values if x == x and x > 0.4]

print(len(all_valid_values))
print(len(validgeqcutoff))
print(len(all_invalid_values))
print(len(invalidgeqcutoff))
"""

sns.boxplot(data=master_df, ax=ax, fliersize=10, showfliers=False, palette=['#9D2727', '#CFB997'])
#sns.violinplot(data=master_df, ax=ax, fliersize=10, showfliers=False, palette=nice_palette)
sns.stripplot(data=master_df, color="black", s=4, ax=ax, jitter=0.2, edgecolor="black", alpha=1, zorder=10, linewidth=0.15)

x=0
add_stat_annotation(ax, data=master_df,
    box_pairs=[(col_names[0],col_names[1])],
    test='Mann-Whitney', text_format='full', loc='outside', verbose=2, text_offset=0.01, line_offset_to_box=(x)/20, line_offset=0.05)

ax.set_facecolor("white")
ax.grid(False)
ax.tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4)
ax.tick_params(axis='y', which='major', left=True, labelsize=10, size=4)

ax.spines['bottom'].set_color('0')
ax.spines['left'].set_color('0')
ax.spines['top'].set_color('1')
ax.spines['right'].set_color('1')
ax.set_ylabel("Methylation at Probe (β)")

for label in ax.get_xticklabels():
    label.set_ha("center")
    label.set_rotation(0)

ax.set_ylim([-0.01, 1.01])
#ax.axhline(0.4, ls="--", lw=0.5, c="black")

fig.tight_layout()

dpi_set = 72
plt.tick_params(bottom='on', left='on')

fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/methylation_new/methylation_by_MAF_plot_reactivated_tRCCs_with_XatRCCs.pdf", dpi=dpi_set, bbox_inches = 'tight')    


