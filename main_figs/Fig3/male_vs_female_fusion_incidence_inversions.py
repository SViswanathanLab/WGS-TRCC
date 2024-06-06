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
sns.set(rc={'figure.figsize':(10,2.2)})
sns.set(font_scale=1.2)
plt.rcParams.update({'font.size': 11})

cs=["#008BB5", "#DACAFF"]

df = pd.read_excel("/Volumes/sviswanathan/projects/TRCC_cohort_summary/all_samples_exomes.xlsx")
df = df[df['Sex']!="Unknown"]
df = df[df['Fusion']!="Unknown"]
df = df[df['CNI (del)'].notna()]

male_full_df = df[df['Sex']=="MALE"]
female_full_df = df[df['Sex']=="FEMALE"]
male_norm = len(male_full_df.index.tolist())
female_norm = len(female_full_df.index.tolist())

print(male_norm)
print(female_norm)

#Partners with at least 5 samples
uq_partner = ['NONO-TFE3', 'RBM10-TFE3']
chr_name = ['chrX', 'chrX']

#uq_partner = ['SFPQ-TFE3','ASPSCR1-TFE3','MED15-TFE3','PRCC-TFE3','LUC7L3-TFE3','Other Autosomal','NONO-TFE3','RBM10-TFE3']
#chr_name = ['chr1', 'chr17', 'chr22', 'chr1', 'chr17', '', 'chrX', 'chrX']

male_count = []
female_count = []
female_cni_del = []
male_cni_del = []
male_cni_gain = []
female_cni_gain = []

for a in uq_partner:
    if a != "Other Autosomal" and a != "Total":
        temp_df = df[df['Fusion']==a]
    elif a == "Total":
        temp_df = df
    else:
        temp_df = df[~(df['Fusion'].isin(uq_partner))]
    male_df = temp_df[temp_df['Sex']=="MALE"]
    female_df = temp_df[temp_df['Sex']=="FEMALE"]
    male_count.append(len(male_df.index.tolist()))
    female_count.append(len(female_df.index.tolist()))

    cni_del_male = male_df[male_df['CNI (del)']=="Yes"]
    cni_gain_male = male_df[male_df['CNI (gain)']=="Yes"]
    cni_del_female = female_df[female_df['CNI (del)']=="Yes"]
    cni_gain_female = female_df[female_df['CNI (gain)']=="Yes"]

    male_cni_del.append(len(cni_del_male.index.tolist()))
    male_cni_gain.append(len(cni_gain_male.index.tolist()))
    female_cni_del.append(len(cni_del_female.index.tolist()))
    female_cni_gain.append(len(cni_gain_female.index.tolist()))

fus_name = []
q=0
while q<len(uq_partner):
    if uq_partner[q] != "Other Autosomal" and uq_partner[q] != "Total":
        fus_name.append(uq_partner[q].split("-")[0] + " (" + chr_name[q] + ", N=" + str(male_count[q]+female_count[q]) + ")")
    else:
        fus_name.append(uq_partner[q].split("-")[0] + " (N=" + str(male_count[q]+female_count[q]) + ")")
    q=q+1

fig, ax = plt.subplots(1, 1)

bar_width = 0.8
x_pos = np.arange(len(uq_partner))

print(female_cni_del)
print(male_cni_del)
print(female_cni_gain)
print(male_cni_gain)

ax.barh(x_pos, -1*np.array(female_count), height=bar_width, color="gainsboro", label='Female', lw=0)
ax.barh(x_pos, np.array(male_count)+5, height=bar_width, color="gainsboro", label='Male', lw=0)
ax.barh(x_pos, -1*np.array(female_cni_gain), height=bar_width, color="red", label='CNI (Gain)', lw=0)
ax.barh(x_pos, 1*np.array(male_cni_gain), height=bar_width, color="red", lw=0)
ax.barh(x_pos, -1*np.array(female_cni_del), left=np.array(female_cni_gain)*-1, height=bar_width, color="blue", label='CNI (Del)', lw=0)
ax.barh(x_pos, 1*np.array(male_cni_del), left=np.array(male_cni_gain), height=bar_width, color="blue", lw=0)

ax.axvline(2.5, lw=85, ls="-", c="white")

#ax.barh(x_pos, -1*np.array(bal_loss), height=bar_width, left=(-1*np.array(unbal)), label='Balanced\nwith Loss', color="#57427A", lw=0)

for i in range(len(x_pos)):
    ev_line_x = male_count[i]*-1
    if i != 6 and i != 7:
        ev_line_x = ev_line_x*2
    print(i)
    print("Value: " + str(ev_line_x))
    #ax.axvline(ev_line_x, ymin=1-i/8-0.035, ymax=7/8-i/8+0.03, color='red', lw=2)

# Set labels and title
#ax.set_ylabel('Fusion Partner',fontsize=19, labelpad=10)
#ax.set_xlabel('Number of TFE3 tRCCs',fontsize=19, labelpad=10)
ax.set_yticks(x_pos)
ax.set_yticklabels(fus_name)
if "Total" in uq_partner:
    ax.set_xticks([-75, -50, -25, 0, 25, 50, 75])
    ax.set_xticklabels([75, 50, 25, 0, 0, 25, 50])
else:
    ax.set_xticks([-15, -10, -5, 0, 5, 10, 15])
    ax.set_xticklabels([15, 10, 5, 0, 0, 5, 10])
#ax.legend(facecolor="white", edgecolor="black", fontsize=19)

ax.set_facecolor("white")
ax.grid(False)
ax.tick_params(axis='x', which='major', bottom=True, labelsize=19, size=4)
ax.tick_params(axis='y', which='major', left=True, labelsize=19, size=4)
ax.spines['bottom'].set_color('0')
ax.spines['left'].set_color('0')
ax.spines['top'].set_color('1')
ax.spines['right'].set_color('1')
ax.invert_yaxis()

fig.tight_layout()
dpi_set = 1200
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/male_vs_female_fusion_partner_incidence_inversions_only.pdf", dpi=dpi_set)
