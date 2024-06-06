import pandas as pd
import collections
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from statannot import add_stat_annotation
import statistics
import scipy
import scipy.stats as stats
from upsetplot import from_indicators, from_memberships, UpSet, plot

sns.set(rc={'figure.figsize':(7,5)})
sns.set(font_scale=1.3)
sns.set(font="Arial")
plt.rcParams.update({'font.size': 12})

df = pd.read_excel("/Volumes/sviswanathan/projects/TRCC_cohort_summary/all_samples_exomes.xlsx")
df = pd.read_excel("/Volumes/sviswanathan/projects/TRCC_cohort_summary/all_samples_exomes.xlsx")
wes_vals_unp = df['WES'].tolist()
fusion_vals_unp = df['Fusion'].tolist()
sex_vals_unp = df['Sex'].tolist()
RNA_vals_unp = df['RNA seq'].tolist()
sample_ids = df['ID'].tolist()

fusion_vals = []
q=0
for a in fusion_vals_unp:
    if a != "Unknown":
        fusion_vals.append(sample_ids[q])
    q=q+1

sex_vals = []
q=0
for a in sex_vals_unp:
    if a != "Unknown":
        sex_vals.append(sample_ids[q])
    q=q+1

wes_vals = []
q=0
for a in wes_vals_unp:
    if a != "NO":
        wes_vals.append(sample_ids[q])
    q=q+1

RNA_vals = []
q=0
for a in RNA_vals_unp:
    if a != "NO":
        RNA_vals.append(sample_ids[q])
    q=q+1

df1_DGEs = wes_vals
df2_DGEs = fusion_vals
df3_DGEs = sex_vals
df4_DGEs = []

dge_1111_unp = list(set(df1_DGEs) & set(df2_DGEs) & set(df3_DGEs) & set(df4_DGEs))
dge_1110_unp = list(set(df1_DGEs) & set(df2_DGEs) & set(df3_DGEs))
dge_1101_unp = list(set(df1_DGEs) & set(df2_DGEs) & set(df4_DGEs))
dge_1011_unp = list(set(df1_DGEs) & set(df3_DGEs) & set(df4_DGEs))
dge_0111_unp = list(set(df2_DGEs) & set(df3_DGEs) & set(df4_DGEs))
dge_1100_unp = list(set(df1_DGEs) & set(df2_DGEs))
dge_1010_unp = list(set(df1_DGEs) & set(df3_DGEs))
dge_1001_unp = list(set(df1_DGEs) & set(df4_DGEs))
dge_0110_unp = list(set(df2_DGEs) & set(df3_DGEs))
dge_0101_unp = list(set(df2_DGEs) & set(df4_DGEs))
dge_0011_unp = list(set(df3_DGEs) & set(df4_DGEs))
dge_1000_unp = df1_DGEs
dge_0100_unp = df2_DGEs
dge_0010_unp = df3_DGEs
dge_0001_unp = df4_DGEs

dge_1111 = dge_1111_unp
dge_1110 = np.setdiff1d(dge_1110_unp,dge_1111_unp)
dge_1101 = np.setdiff1d(dge_1101_unp,dge_1111_unp)
dge_1011 = np.setdiff1d(dge_1011_unp,dge_1111_unp)
dge_0111 = np.setdiff1d(dge_0111_unp,dge_1111_unp)
dge_1100 = np.setdiff1d(dge_1100_unp,dge_1110_unp+dge_1101_unp)
dge_1010 = np.setdiff1d(dge_1010_unp,dge_1110_unp+dge_1011_unp)
dge_1001 = np.setdiff1d(dge_1001_unp,dge_1101_unp+dge_1011_unp)
dge_0110 = np.setdiff1d(dge_0110_unp,dge_0111_unp+dge_1110_unp)
dge_0101 = np.setdiff1d(dge_0101_unp,dge_0111_unp+dge_1101_unp)
dge_0011 = np.setdiff1d(dge_0011_unp,dge_0111_unp+dge_1011_unp)
dge_1000 = np.setdiff1d(dge_1000_unp,dge_0100_unp+dge_0010_unp+dge_0001_unp)
dge_0100 = np.setdiff1d(dge_0100_unp,dge_1000_unp+dge_0010_unp+dge_0001_unp)
dge_0010 = np.setdiff1d(dge_0010_unp,dge_0100_unp+dge_1000_unp+dge_0001_unp)
dge_0001 = np.setdiff1d(dge_0001_unp,dge_0100_unp+dge_0010_unp+dge_1000_unp)
#['1111_FUUR1_UOK146_UOK109_sTFE', '1110', '1101', '1011', '0111', '1100', '1010', '1001', '0110', '0101', '0011', '1000', '0100', '0010', '0001']

len_vals = [len(dge_1111), len(dge_1110), len(dge_1101), len(dge_1011), len(dge_0111), len(dge_1100), len(dge_1010), len(dge_1001), len(dge_0110), len(dge_0101), len(dge_0011), len(dge_1000), 0, len(dge_0010), len(dge_0001)]

membership_list = [
['WES/WGS','Fusion Call','Sex Annotation','RNA-seq'],
['WES/WGS','Fusion Call','Sex Annotation'],
['WES/WGS','Fusion Call','RNA-seq'],
['WES/WGS','Sex Annotation','RNA-seq'],
['Fusion Call','Sex Annotation','RNA-seq'],
['WES/WGS','Fusion Call'],
['WES/WGS','Sex Annotation'],
['WES/WGS','RNA-seq'],
['Fusion Call','Sex Annotation'],
['Fusion Call','RNA-seq'],
['Sex Annotation','RNA-seq'],
['WES/WGS'],
['Fusion Call'],
['Sex Annotation'],
['RNA-seq'],
]

# Create two new lists removing values if their corresponding length is 0
new_len_vals, new_membership_list = zip(*[(l, m) for l, m in zip(len_vals, membership_list) if l > 0])

# Print the new lists
print(new_len_vals)
print(new_membership_list)

example = from_memberships(new_membership_list, data=new_len_vals)

plot(example)
ax1 = plt.gca()

ax1.set_xlabel('Number of TFE3 tRCCs', fontsize=13)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.xaxis.set_tick_params(labelsize=10)
ax1.set_facecolor("white")
ax1.grid(False)
ax1.tick_params(axis='y', which='major', left=True, labelleft=True, labelsize=13, size=4)
ax1.tick_params(axis='x', which='major', bottom=True, labelbottom=True, labelsize=13, size=4)
ax1.spines['bottom'].set_color('0')
ax1.spines['left'].set_color('0')

for label in ax1.get_xticklabels():
    label.set_ha("center")
    label.set_rotation(0)

ax1.set_ylim([0, 150])

plt.tight_layout()
dpi_set = 1200

plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/cohort_composition_upset.pdf", dpi=dpi_set, bbox_inches = 'tight')


"""



upset_data = pd.DataFrame(data)
print(upset_data)

#fig, ax = plt.subplots(figsize=(8, 6))
#UpSet(matrix, subset_size="sum")
#plt.show()
"""
"""
#ax1.barh(uq_cohorts, samp1, color='black', lw=0)
#ax1.barh(uq_cohorts, samp2, left=np.array(samp1), color='red', lw=0)
ax1.set_xlabel('Number of TFE3 tRCCs', fontsize=13)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.xaxis.set_tick_params(labelsize=10)
ax1.set_facecolor("white")
ax1.grid(False)
ax1.tick_params(axis='y', which='major', left=True, labelleft=True, labelsize=13, size=4)
ax1.tick_params(axis='x', which='major', bottom=True, labelbottom=True, labelsize=13, size=4)
ax1.spines['bottom'].set_color('0')
ax1.spines['left'].set_color('0')

for label in ax1.get_xticklabels():
    label.set_ha("center")
    label.set_rotation(0)

fig.tight_layout()
dpi_set = 1200

plt.subplots_adjust(hspace=0)

plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/cohort_composition_upset.pdf", dpi=dpi_set, bbox_inches = 'tight')
"""