pip install comut

## import requirements
from comut import comut
from comut import fileparsers
import pandas as pd
import numpy as np
import palettable
%matplotlib inline
%config InlineBackend.figure_format = 'retina'
import matplotlib
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator # this function sets the location of the minor tick mark
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import date

sample_sorted_df = pd.read_csv ('sample_sorted_rev1.csv', sep = ",") # sorted by sample wise
purity_df = pd.read_csv('purity.csv', sep = ',')
ploidy_df = pd.read_csv('ploidy_manual.csv', sep = ',')
## for one or more samples
tum_df = pd.read_csv('tum_data_rev1.csv', sep =",")
tmb_df = pd.read_csv('tmb_merg_rev1.csv', sep=",")
# add fusion partner info
fusion_df = pd.read_csv ('fusion_partner.csv', sep = ",")
## updated genes for revision 1
genes = ['UBR5', 'PTPRK', 'AR', 'IL21R', 'LEF1', 'MAML2', 'CDCA2', 'CEP135', 'SETX', 'TYRO3', 'WDR36', 'XPC', 'IKZF3', 'TP53', 'HNF1A', 'KDM5C', 'SETD2', 'KMT2C','CDKN2C', 'BRCA1', 'ATM', 'ARID1A', 'NONE']

# set sample order
#sample_order = list(sample_sorted_df['sample'])
## add mapping and kwargs for colors, range for the value data range
## color mapping for categorical and continuous data

purity_mapping = {'Purity': 'rosybrown'}
ploidy_mapping = {'2N': 'slateblue', '4N': 'lightsalmon'}
sample_sorted_mapping = {'Frozen-WGS':'#2A385B','Frozen-LR':'#DECBE3', 'FFPE-WGS':'#8B7F47'}

indicator_kwargs = {'color': 'black', 'marker': 'o', 'linewidth': 1, 'markersize': 3}

fusion_mapping = {'X-autosome':'cadetblue', 'X-inv':'goldenrod'}
tum_mapping ={'None': 'white','Missense': 'yellow', 'Frameshift': 'purple', 'Splice site': 'darkseagreen', 'Stop gained': 'magenta'}
tmb_bar_mapping = {'INDEL':'darkorange','SNV':'teal'}
value_order = ['Frameshift']
priority = ['Frameshift']

wgs_comut = comut.CoMut()
wgs_comut.add_categorical_data(tum_df, name = 'Mutation type', category_order = genes, value_order= value_order, mapping = tum_mapping, tick_style = 'italic', priority=priority)
wgs_comut.add_categorical_data(sample_sorted_df, name = 'Sample type', mapping = sample_sorted_mapping)
wgs_comut.add_categorical_data(fusion_df, name = 'Fusion partner', mapping = fusion_mapping)
wgs_comut.add_categorical_data(ploidy_df, name = 'Ploidy', mapping = ploidy_mapping)
wgs_comut.add_bar_data(purity_df, name = 'Purity', mapping = purity_mapping, stacked = True, bar_kwargs = {'width': 0.8, 'edgecolor': 'none'}, ylabel = 'Purity')
wgs_comut.add_bar_data(tmb_df,   name = 'Mutation rate', mapping = tmb_bar_mapping, stacked = True, ylabel = 'Muts/Mb', bar_kwargs = {'width': 0.8, 'edgecolor': 'none'})
custom_rcParams = {'font.family': 'arial', 'font.size': 4.7}
rcParams.update(custom_rcParams)
border_white = ['border_white']

wgs_comut.plot_comut(x_padding = 0.04, y_padding = 0.05, wspace=0.05)
heights = {'Mutation rate':7.5, 'Purity':3.5, 'Ploidy':1.2, 'Sample type':1.2}
border_white = ['border_white']
x_padding = 0.04 # the x distance between patches in comut
y_padding = 0.04 # the y distance between patches in comut
heights = {'Mutation rate':7.5, 'Purity':3.5, 'Ploidy':1.2, 'Sample type':1.2}
wgs_comut.add_unified_legend(border_white = border_white, ncol = 1)
wgs_comut.axes['Mutation rate'].set_ylim(0,5)
wgs_comut.axes['Mutation rate'].set_ylabel('Mutations/Mb', rotation = 'horizontal', horizontalalignment='right', x=-9)
wgs_comut.axes['Purity'].set_ylim(0,1)
wgs_comut.axes['Purity'].set_ylabel('Purity', rotation = 'horizontal', horizontalalignment='right', x=-4)
fig = plt.gcf()

#wgs_comut.add_axis_legend(border_white = border_white, ncol = 3)
#update fonts
custom_rcParams = {'font.family': 'arial',
                   'font.size': 4.7, 'legend.fontsize' : 4.7}
#rcParams.update(custom_rcParams)
border_white = ['border_white']
x_padding = 0.04 # the x distance between patches in comut
y_padding = 0.04 # the y distance between patches in comut
#plot
wgs_comut.plot_comut(x_padding = x_padding, y_padding = y_padding, figsize=(8, 3), wspace = 0.12, hspace = 0.4)
wgs_comut.axes['Mutation rate'].set_ylim(0,6)
wgs_comut.axes['Mutation rate'].set_yticks([0,5])
wgs_comut.axes['Mutation rate'].set_ylabel('Mutations/Mb', rotation = 'horizontal', horizontalalignment='right', x=-9)
wgs_comut.axes['Purity'].set_ylim(0,1)
wgs_comut.axes['Purity'].set_yticks([0,0.5,1])
wgs_comut.axes['Purity'].set_ylabel('Purity', rotation = 'horizontal', horizontalalignment='right', x=-4)

fig = plt.gcf()

fig.set_size_inches(8,3)
wgs_comut.add_unified_legend(border_white = border_white, ncol = 1)
