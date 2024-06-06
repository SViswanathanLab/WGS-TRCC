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
sns.set(rc={'figure.figsize':(4,6)})
sns.set(font_scale=1)
plt.rcParams.update({'font.size': 13})

cs=['#008BB5', '#DACAFF']

df = pd.read_excel("/Volumes/sviswanathan/projects/TRCC_cohort_summary/all_samples_exomes.xlsx")
df = df[df['Sex']!="Unknown"]

#df = df[df['Unbalanced'].notna()]
#df = df[df['Fusion']=="SFPQ-TFE3"]

male_df = df[df['Sex']=="MALE"]
female_df = df[df['Sex']=="FEMALE"]

# create a figure with one subplot
fig, ax = plt.subplots(1, 1)

#Plot pie chart here:

dpi_set = 1200

recipe = ["Balanced",
          "Unbalanced"]

data = [len(male_df.index.tolist()), len(female_df.index.tolist())]
print(data[0])
print(data[1])

#Out of known samples

wedges, texts = ax.pie(data, wedgeprops=dict(width=0.5), startangle=21, pctdistance=0.74, textprops={'fontsize': 5, 'rotation':0}, colors=cs)

bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"),
          bbox=bbox_props, zorder=3, va="center")

for i, p in enumerate(wedges):
    
    if i != 2:
   
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        #ax.annotate(recipe[i], xy=(x, y), xytext=(1.1*np.sign(x), 1.2*y),
        #            horizontalalignment=horizontalalignment, **kw)
    else:
   
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(48)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        #ax.annotate(recipe[i], xy=(x, y), xytext=(1.1*np.sign(x), 1.2*y+0.4),
        #            horizontalalignment=horizontalalignment, **kw)

"""
for patch, txt in zip(wedges, percent):
    if (patch.theta2 - patch.theta1) <= 30:
        # the angle at which the text is normally located
        angle = (patch.theta2 + patch.theta1) / 2.
        # new distance to the pie center
        #x = patch.r * 1.2 * np.cos(angle * np.pi / 180)
        #y = patch.r * 1.2 * np.sin(angle * np.pi / 180)
        # move text to new position
        #txt.set_position((1.07*x, 0.92*y))
        #txt.set_color('black')
"""

ax.set_facecolor("white")
ax.grid(False)
ax.tick_params(axis='x', which='major', bottom=True, labelsize=13, size=4)
ax.tick_params(axis='y', which='major', left=True, labelsize=13, size=4)
ax.spines['bottom'].set_color('0')
ax.spines['left'].set_color('0')
ax.spines['top'].set_color('1')
ax.spines['right'].set_color('1')

fig.tight_layout()
dpi_set = 1200
fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/male_vs_female_donut.pdf", dpi=dpi_set)
