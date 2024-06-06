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

sns.set(rc={'figure.figsize':(3.75*1.4+1,3.6*1.4)})
sns.set(font_scale=1.2)
sns.set(font="Arial")
plt.rcParams.update({'font.size': 18})

df = pd.read_excel("/Volumes/sviswanathan/projects/TRCC_cohort_summary/all_samples_exomes.xlsx")
df = df[df['Partner_chromosome']!="Unknown"]
df = df[df['Partner_chromosome']!="Undetected"]

uq_cohorts = sorted(list(set(df['Dataset'].tolist())))

# get the total number of samples in each cohort
total_samples_count = []
for cohort_uq in uq_cohorts:
    temp_df = df[df['Dataset']==cohort_uq]
    total_samples_count.append(len(temp_df.index.tolist()))

new_cohorts = [x for _, x in sorted(zip(total_samples_count, uq_cohorts))]
total_samples_count = sorted(total_samples_count)
uq_cohorts = new_cohorts
total_samples_count = total_samples_count

# plot the total number of samples in each cohort as a bar chart in the upper subplot
fig, ax2 = plt.subplots(1, 1)

"""
ax1.barh(uq_cohorts, total_samples_count, color='black')
ax1.set_xlabel('Number of tRCCs', fontsize=13)
ax1.spines['bottom'].set_visible(True)
ax1.spines['left'].set_visible(False)
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
"""

uq_cohorts = sorted(list(set(df['Dataset'].tolist())))
total_male_chrX_count = []
total_female_chrX_count = []
total_male_auto_count = []
total_female_auto_count = []

for cohort_uq in uq_cohorts:

    temp_df = df[df['Dataset']==cohort_uq]

    male_df = temp_df[temp_df['Sex']=="MALE"]
    female_df = temp_df[temp_df['Sex']=="FEMALE"]

    male_chrX = male_df[male_df['Partner_chromosome']=="chrX"]
    male_autosomal = male_df[male_df['Partner_chromosome']!="chrX"]

    female_chrX = female_df[female_df['Partner_chromosome']=="chrX"]
    female_autosomal = female_df[female_df['Partner_chromosome']!="chrX"]

    male_chrX_count = len(male_chrX.index.tolist())
    male_autosomal_count = len(male_autosomal.index.tolist())
    female_chrX_count = len(female_chrX.index.tolist())
    female_autosomal_count = len(female_autosomal.index.tolist())

    total_male_auto_count.append(male_autosomal_count)
    total_female_auto_count.append(female_autosomal_count)
    total_male_chrX_count.append(male_chrX_count)
    total_female_chrX_count.append(female_chrX_count)

print(sum(total_male_auto_count))
print(sum(total_female_auto_count))
print(sum(total_male_chrX_count))
print(sum(total_female_chrX_count))

output_df = pd.DataFrame([total_male_chrX_count,total_female_chrX_count,total_male_auto_count,total_female_auto_count]).T
output_df.columns = ['Xinv-M','Xinv-F','X:A-M','X:A-F']
output_df.index = uq_cohorts
output_df.loc["Total"] = output_df.sum()

values = output_df.values.tolist()
zeros_list = []
totals = []
for a in values:
    totals.append(sum(a))
    zeros_list.append(0)

output_df['Total'] = totals

master_df = output_df
# Normalize cols 1-2 separately
cols12_norm = master_df.iloc[:, :2].div(master_df.iloc[:, :2].sum(axis=1), axis=0)
# Normalize cols 3-4 separately
cols34_norm = master_df.iloc[:, 2:4].div(master_df.iloc[:, 2:4].sum(axis=1), axis=0)

master_df['Zero'] = zeros_list
col5_norm = master_df['Zero']
del master_df['Zero']

# Concatenate the normalized dataframes and the unnormalized columns
master_df = pd.concat([cols12_norm, cols34_norm, col5_norm, master_df.iloc[:, 5:]], axis=1)

ax2.set_facecolor("white")
ax2.grid(False)
ax2.tick_params(axis='x', which='major', top=False, labeltop=False, bottom=True, labelbottom=True, labelsize=18, size=4)
ax2.tick_params(axis='y', which='major', left=True, labelsize=16, size=4)
ax2.spines['bottom'].set_color('0')
ax2.spines['left'].set_color('0')

width = 0.38

#ax.set_xlabel("hg38 chrX Position", labelpad=10)
ax2.set_ylabel("Number of $\it{TFE3}$ tRCCs",fontsize=16)
ind = np.arange(0, len(uq_cohorts))
ax2.set_xticks(ind)
ax2.set_xticklabels(uq_cohorts[::-1], fontsize=18)

col_names = ['chrX\nInversion\nMale','chrX\nInversion\nFemale','X:Autosome\nTranslocation\nMale','X:Autosome\nTranslocation\nFemale','Total']
output_df.columns = col_names
master_df.columns = col_names

nice_palette = ['#47B5DC', '#F2E4FF', '#00709C', '#BFA3D9', '#DDDD00', '#36E2EE', '#008BB5', '#DACAFF', '#FFCAFF']
nice_palette = ['#008BB5', '#DACAFF','#008BB5', '#DACAFF']

master_df = master_df.iloc[[len(master_df.index.tolist())-1]]
output_df = output_df.iloc[[len(output_df.index.tolist())-1]]
master_df.index = ['Aggregate Cohort']
output_df.index = ['Aggregate Cohort']

a,p = scipy.stats.chisquare([21, 8], [29*1/2,29*1/2])
print(p)

a,p = scipy.stats.chisquare([51, 102], [153*1/2,153*1/2])
print(p)

a,p = scipy.stats.chisquare([11, 7], [9,9])
print(p)



ax2.bar(output_df.columns.tolist()[0:4][::-1], output_df.values.tolist()[0][0:4][::-1], width=0.9, lw=0, color=nice_palette[0:4][::-1])
ax2.set_xticks([0, 1, 2, 3])
ax2.set_xticklabels(col_names[0:4])

fig.tight_layout()
dpi_set = 1200

plt.subplots_adjust(hspace=0)

plt.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/tRCC_male_female_inversion_freq_cohort_level.pdf", dpi=dpi_set, bbox_inches = 'tight')
