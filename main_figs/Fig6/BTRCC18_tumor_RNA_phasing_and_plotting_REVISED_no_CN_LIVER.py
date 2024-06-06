import glob, os
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
from os import listdir
import gzip
from statsmodels.stats.proportion import proportions_ztest
from matplotlib.collections import BrokenBarHCollection

plt.rcParams["font.family"] = "Arial"
sns.set(rc={'figure.figsize':(11.5,11.75)})
sns.set(font_scale=1)
plt.rcParams.update({'font.size': 12})

hapbased_phase = pd.read_csv("/Users/ananthansadagopan/Downloads/finalHap.txt", sep="\t")

reverse_val = False

def split_advanced(strng, sep, pos):
    strng = strng.split(sep)
    return sep.join(strng[:pos]), sep.join(strng[pos:])

def get_vcf_names(vcf_path):
    header = []
    with gzip.open(vcf_path, "rt") as ifile:
        for line in ifile:
            if line.startswith("##"):
                header.append(line)
            if line.startswith("#CHROM"):
                vcf_names = [x.split("\n")[0] for x in line.split('\t')]
                break
    ifile.close()
    return vcf_names, header

dpi_set = 1200



"""
STEP 1:

PHASING RNA
"""




#Statistical phasing calls on common germline variants left of TFE3
df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/Xi_Xa_fusion/BTRCC18_phased/BTRCC18_phased_germline_vcf.csv")
df = df[df['POS']<49043410]
pos_list = df['POS'].tolist()
phase_list = df['Phase'].tolist()

#File used to correct switching errors on the above calls
df2 = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/Xi_Xa_fusion/BTRCC18_phased/hap_fraction.xlsx")

hap = df2['hap'].tolist()
hap_dict = dict(zip(df2.index.tolist(), hap))
new_phase = []

a=0
while a<len(pos_list):
    curr_pos = math.floor(pos_list[a]/10000)
    curr_phase = phase_list[a]
    curr_corr = hap_dict[curr_pos]
    if int(curr_corr) == -1:
        new_phase.append(curr_phase)
    elif curr_phase == "1|0":
        new_phase.append("0|1")
    elif curr_phase == "0|1":
        new_phase.append("1|0")
    else:
        new_phase.append("Unknown")
    a=a+1

df['new_phase'] = new_phase
variant_pos = df['POS'].tolist()

#Now use allelic depth calls to phase the right of TFE3
df3 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/Xi_Xa_fusion/BTRCC18_phased/Xq_phase_from_allelic_depth_DNA_revised_clonal_beforeWGD.csv")
df3 = df3[df3['POS']>=49043410]
new_pos = df3['POS'].tolist()
phase2 = df3['Phase'].tolist()

#Combine the two sets of phasing calls together
pos_dict = dict(zip(variant_pos+new_pos, new_phase+phase2))

#Assign phases to the RNA ASE output
df3 = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/RNAseq_paired_all_chr/allchr_RNA_paired_on_germline_hets_unphased_LiverTum_2.log", sep="\t")
df3 = df3[df3['contig']=="chrX"]
pos_list = df3['position'].tolist()
final_phases = []

for a in pos_list:
    try:
        curr_phase = pos_dict[a]
    except:
        curr_phase = "Unknown"
    final_phases.append(curr_phase)

df3['Phase'] = final_phases
df2 = df3[df3['Phase']!="Unknown"]




"""
STEP 2:

DETERMINING DNA COUNTS FOR PHASED RNA HETS
"""



df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/Xi_Xa_fusion/BTRCC18_phased/merged_chrX_variants_common_rare_somatic.csv")

positions = df2['position'].tolist()
phase_vals = df2['Phase'].tolist()

pos_vals = df['POS'].tolist()
vals = df['B_TRCC_18_Tumor'].tolist()

new_vals1 = []
new_vals2 = []
new_vals3 = []

for a in vals:
    new_vals3.append(a.split(":")[0])
    sub = a.split(":")[1].split(",")
    if sub != ['.']:
        new_vals1.append(sub[0])
        new_vals2.append(sub[1])
    else:
        new_vals1.append(".")
        new_vals2.append(".")

dict1 = dict(zip(pos_vals, new_vals1))
dict2 = dict(zip(pos_vals, new_vals2))
dict3 = dict(zip(pos_vals, new_vals3))

vals1 = []
vals2 = []
call = []

a=0
while a<len(positions):
    call.append(dict3[int(positions[a])])
    vals1.append(dict1[int(positions[a])])
    vals2.append(dict2[int(positions[a])])
    a=a+1

df2['DNA_ref_count'] = vals1
df2['DNA_alt_count'] = vals2
df2['call'] = call

df2 = df2.replace('.',0)
df2['total_DNA_Count'] = df2['DNA_ref_count'].astype(int) + df2['DNA_alt_count'].astype(int)
#df2 = df2[df2['call']!="0/0"]

vals = df['LungPlTum_1'].tolist()

new_vals1 = []
new_vals2 = []
new_vals3 = []

for a in vals:
    new_vals3.append(a.split(":")[0])
    sub = a.split(":")[1].split(",")
    if sub != ['.']:
        new_vals1.append(sub[0])
        new_vals2.append(sub[1])
    else:
        new_vals1.append(".")
        new_vals2.append(".")

dict1 = dict(zip(pos_vals, new_vals1))
dict2 = dict(zip(pos_vals, new_vals2))
dict3 = dict(zip(pos_vals, new_vals3))

vals1 = []
vals2 = []
call = []

a=0
while a<len(positions):
    call.append(dict3[int(positions[a])])
    vals1.append(dict1[int(positions[a])])
    vals2.append(dict2[int(positions[a])])
    a=a+1

df2['DNA_ref_count_Lung'] = vals1
df2['DNA_alt_count_Lung'] = vals2
df2['call_Lung'] = call

df2.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/Xi_Xa_fusion/BTRCC18_phased/phased_ASE_LiverTum2_chrX_revised_RNA_and_DNA_0310.csv", index=False)

input_counts_with_phase_df = df2

#Compiled seg file
df_abs_cn = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/BTRCC18_allelic_CN.xlsx")

#Defining escapee exon boundaries (creating list of excluded regions)
chrX_exons_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/chrX_all_exons_gencode_v43_UCSC.csv")
gene_classes_df = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/kidney_specific_chrX_gene_class_list_cotton.xlsx")
escape_gene = [x for x in gene_classes_df['Escape_plus_Variable'].tolist() if x == x]
escape_df = chrX_exons_df[chrX_exons_df['name2'].isin(escape_gene)]
exon_starts = escape_df['exonStarts'].tolist()
exon_ends = escape_df['exonEnds'].tolist()

#Create exclusion lists with PARs already included
new_starts = [1, 155701383]
new_ends = [2781479, 156030895]
for a in exon_starts:
    temp = a.split(",")
    for b in temp:
        new_starts.append(b)
for a in exon_ends:
    temp = a.split(",")
    for b in temp:
        new_ends.append(b)

valid_new_starts = []
valid_new_ends = []
q=0
while q<len(new_starts):
    if new_starts[q] == "" or new_ends[q] == "":
        pass
    else:
        valid_new_starts.append(new_starts[q])
        valid_new_ends.append(new_ends[q])
    q=q+1

KS_pval = []
vals = ['BTRCC18']

#Iterate over each file
for www in vals:
    value = "B_TRCC_18_Tumor"
    temp_chr = "chrX"
    temp_id = '18'

    #Find germline hets with DNA VAF within 1SD of the mean across chrX
    sd_cutoff = 1
    print("reading VCF")

    all_files_in_dir = os.listdir("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/germline_hets/")
    temp_DNA_id = [x for x in all_files_in_dir if temp_id in x]
    print(temp_DNA_id)

    path1 = "/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/germline_hets/" + temp_DNA_id[0]
    names, header_temp = get_vcf_names(path1)
    df = pd.read_csv(path1, comment='#', delim_whitespace=True, header=None, names=names, compression="gzip")
    df = df[df['#CHROM']==temp_chr]
    print("read VCF")
    
    normals = [temp_id]
    ref_df_rows = []
    alt_df_rows = []

    for a in normals:
        temp_val = df["CELL"].tolist()
        new_val1 = []
        new_val2 = []
        for q in temp_val:
            temp = q.split(":")[1]
            temp2 = q.split(":")[0]
            sub = temp.split(",")
            if sub != ['.'] and temp2 != "0/0" and temp2 != "1/1":
                new_val1.append(sub[0])
                new_val2.append(sub[1])
            else:
                new_val1.append(float("nan"))
                new_val2.append(float("nan"))          
        ref_df_rows.append(new_val1)
        alt_df_rows.append(new_val2)

    print("converting rows")
    ref_df = pd.DataFrame(ref_df_rows).T
    alt_df = pd.DataFrame(alt_df_rows).T
    ref_df_rows = ref_df.values.tolist()
    alt_df_rows = alt_df.values.tolist()

    AF = []
    total_sum = []

    q=0
    while q<len(ref_df_rows):
        temp_row1 = [int(x) for x in ref_df_rows[q] if x == x]
        temp_row2 = [int(x) for x in alt_df_rows[q] if x == x]
        try:
            sum1 = sum(temp_row1)/len(temp_row1)
        except:
            sum1 = float("nan")
        try:
            sum2 = sum(temp_row2)/len(temp_row2)
        except:
            sum2 = float("nan")
        if sum1 == 0 and sum2 == 0:
            AF.append(float("nan"))
            total_sum.append(0)
        else:
            AF.append(sum1/(sum1 + sum2))
            total_sum.append(sum1+sum2)
        q=q+1
        
    df['AF_normal'] = AF
    df['total_sum'] = total_sum
    df = df[df['AF_normal'].notna()]
    df = df[df['total_sum'].notna()]

    #Find germline hets with at least 10 supporting DNA reads
    df = df[df['total_sum']>=10]
    curr_row = df['AF_normal'].tolist()
    variant_pos = df['POS'].tolist()

    mean = sum(curr_row) / len(curr_row) 
    variance = sum([((x - mean) ** 2) for x in curr_row]) / (len(curr_row)-1) 
    standard_dev = variance**0.5

    ys_temp = []
    for t in curr_row:
        z_score = (t-mean)/standard_dev
        ys_temp.append(z_score)

    valid_pos = []

    a=0
    while a<len(variant_pos):
        if ys_temp[a] <= sd_cutoff*-1 or ys_temp[a] >= sd_cutoff:
            pass
        else:
            valid_pos.append(variant_pos[a])
        a=a+1

    xs = []

    #Plot ideogram
    color_lookup = {
                    'gneg': (.7, .7, .7),
                    'gpos25': (.6, .6, .6),
                    'gpos50': (.4, .4, .4),
                    'gpos75': (.2, .2, .2),
                'gpos100': (0., 0., 0.),
                    'acen': (.8, .4, .4),
                    'gvar': (.8, .8, .8),
                    'stalk': (.9, .9, .9),
                }

    height = 1
    spacing = 0

    def ideograms(fn):
        last_chrom = None
        fin = open(fn)
        fin.readline()
        xranges, colors = [(0,4400000)], [(.7, .7, .7)]
        ymin = 0

        for line in fin:
            chrom, start, stop, label, stain = line.strip().split('\t')
            start = int(start)
            stop = int(stop)
            width = stop - start
            if chrom == last_chrom or (last_chrom is None):
                xranges.append((start, width))
                colors.append(color_lookup[stain])
                last_chrom = chrom
                continue

            ymin += height + spacing
            yrange = (ymin, height)
            yield xranges, yrange, colors, last_chrom
            xranges, colors = [], []
            xranges.append((start, width))
            colors.append(color_lookup[stain])
            last_chrom = chrom

        ymin += height + spacing
        yrange = (ymin, height)
        yield xranges, yrange, colors, last_chrom

    fig, axs = plt.subplots(7, 3, sharey='row', gridspec_kw={'height_ratios': [9,9,9,9,6,1,1], 'width_ratios': [1, 3, 1]})

    #axs[4,1] to axs[3,1], axs[5,1] to axs[2,1], axs[3,1] to axs[6,1], axs[2,1] to axs[0,1], axs[6,1] to axs[1,1], axs[0,1] to axs[4,1] and axs[1,1] to axs[5,1]

    d = {}
    yticks = []
    yticklabels = []

    # ideogram.txt downloaded from UCSC's table browser
    for xranges, yrange, colors, label in ideograms('/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/cytoBandIdeo_chrX.txt'):
        print(xranges)
        coll = BrokenBarHCollection(xranges, yrange, facecolors=colors, lw=0)
        axs[5,1].add_collection(coll)
        center = yrange[0] + yrange[1]/2.
        yticks.append(center)
        yticklabels.append(label)
        d[label] = xranges

    axs[5,1].axis('tight')
    axs[5,1].set_yticks([])
    axs[5,1].set_yticklabels([])
    axs[5,1].set_xticks([])

    axs[5,1].set_xlim([-1000000, 157040895])

    # ideogram.txt downloaded from UCSC's table browser
    for xranges, yrange, colors, label in ideograms('/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/cytoBandIdeo_chrX.txt'):
        print(xranges)
        coll = BrokenBarHCollection(xranges, yrange, facecolors=colors, lw=0)
        axs[6,1].add_collection(coll)
        center = yrange[0] + yrange[1]/2.
        yticks.append(center)
        yticklabels.append(label)
        d[label] = xranges

    axs[6,1].axis('tight')
    axs[6,1].set_yticks([])
    axs[6,1].set_yticklabels([])
    axs[6,1].set_xticks([])

    axs[6,1].set_xlim([-1000000, 157040895])

    #Plot copy number seg for Primary
    #axs[3,1].set_xlim([-1000000, 157040895])
    #axs[3,1].set_ylim([-0.1, 2.1]) 
    #axs[3,1].set_yticks([0, 1, 2]) 
    #axs[3,1].set_yticklabels([0, 1, 2]) 

    #Purity: chr9:21700000-22000000
    df_pur = df_abs_cn[df_abs_cn['contig']=="hs9"]
    df_pur = df_pur[df_pur['start']>=21700000]
    df_pur = df_pur[df_pur['start']<=21900000]
    df_pur = df_pur[['LiverTum2_A','LiverTum2_B']]
    vals = df_pur.values.tolist()
    sum_list = []
    for tempX in vals:
        sum_list.append(sum(tempX))
    normal_contam_estimate = sum(sum_list)/len(sum_list)
    #purity = 1 - normal_contam_estimate
    purity = 0.80

    df_af = df_abs_cn[df_abs_cn['contig']=="hsX"]
    df_af = df_af[df_af['start']>=2800000]
    df_af = df_af[df_af['start']<=155600000]
    temp1 = list(df_af['LiverTum2_A'])
    temp2 = list(df_af['LiverTum2_B'])

    ploidy = 1.96

    new_temp1 = []
    new_temp2 = []
    for qq in temp1:
        if qq >= 0:
            new_temp1.append(qq)
            #new_temp1.append((qq*(ploidy*purity+2*(1-purity))-(1-purity))/purity)
        else:
            new_temp1.append(0)
    
    for qq in temp2:
        if qq >= 0:
            new_temp2.append(qq)
            #new_temp2.append((qq*(ploidy*purity+2*(1-purity))-(1-purity))/purity)
        else:
            new_temp2.append(0)

    invalid_interval = []
    meantemp1 = sum(new_temp1)/len(new_temp1)
    stdtemp1 = statistics.stdev(new_temp1)
    meantemp2 = sum(new_temp2)/len(new_temp2)
    stdtemp2 = statistics.stdev(new_temp2)

    start_vals = df_af['start'].tolist()
    end_vals = df_af['end'].tolist()
    invalid_starts = []
    invalid_ends = []
    qx=0
    while qx<len(start_vals):
        tempval=0
        for postemp in xs:
            if postemp >= float(start_vals[qx]) and postemp < float(end_vals[qx]):
                tempval=1
        if tempval == 0:
            #pass
            if start_vals[qx]>49028726:
                axs[1,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp1[qx]), float(new_temp1[qx])), marker=None, linestyle="-", lw=3, c="black")
                axs[1,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp2[qx]), float(new_temp2[qx])), marker=None, linestyle="-", lw=3, c="orange", alpha=0.6)
            else:
                axs[1,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp1[qx]), float(new_temp1[qx])), marker=None, linestyle="-", lw=3, c="lightgray")
                axs[1,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp2[qx]), float(new_temp2[qx])), marker=None, linestyle="-", lw=3, c="lightgray", alpha=0.6)                
        else:
            #pass
            axs[1,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp1[qx]), float(new_temp1[qx])), marker=None, linestyle="-", lw=3, c="green", zorder=10)
            axs[1,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp2[qx]), float(new_temp2[qx])), marker=None, linestyle="-", lw=3, c="green", zorder=10)  
        if new_temp1[qx] < meantemp1-stdtemp1 or new_temp1[qx] > meantemp1+stdtemp1 or new_temp2[qx] < meantemp2-stdtemp2 or new_temp2[qx] > meantemp2+stdtemp2: #Identify sites with CN imbalance
            invalid_starts.append(start_vals[qx])
            invalid_ends.append(end_vals[qx])
        qx=qx+1

    #Plot copy number seg for Lung Pl #1
    axs[2,1].set_xlim([-1000000, 157040895])
    #axs[2,1].set_ylim([-0.1, 2.1]) 
    #axs[2,1].set_yticks([0, 1, 2]) 
    #axs[2,1].set_yticklabels([0, 1, 2]) 

    df_pur = df_abs_cn[df_abs_cn['contig']=="hs9"]
    df_pur = df_pur[df_pur['start']>=21700000]
    df_pur = df_pur[df_pur['start']<=21900000]
    df_pur = df_pur[['LungPl1_A','LungPl1_B']]
    vals = df_pur.values.tolist()
    sum_list = []
    for tempX in vals:
        sum_list.append(sum(tempX))
    normal_contam_estimate = sum(sum_list)/len(sum_list)
    #purity = 1 - normal_contam_estimate
    purity = 0.5446

    df_af = df_abs_cn[df_abs_cn['contig']=="hsX"]
    df_af = df_af[df_af['start']>=2800000]
    df_af = df_af[df_af['start']<=155600000]
    temp1 = list(df_af['LungPl1_A'])
    temp2 = list(df_af['LungPl1_B'])

    ploidy = 1.99

    new_temp1 = []
    new_temp2 = []
    for qq in temp1:
        if qq >= 0:
            new_temp1.append(qq)
            #new_temp1.append((qq*(ploidy*purity+2*(1-purity))-(1-purity))/purity)
        else:
            new_temp1.append(0)
    
    for qq in temp2:
        if qq >= 0:
            new_temp2.append(qq)
            #new_temp2.append((qq*(ploidy*purity+2*(1-purity))-(1-purity))/purity)
        else:
            new_temp2.append(0)

    start_vals = df_af['start'].tolist()
    end_vals = df_af['end'].tolist()
    qx=0
    while qx<len(start_vals):
        if start_vals[qx]>49028726:
            axs[0,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp1[qx]), float(new_temp1[qx])), marker=None, linestyle="-", lw=3, c="black")
            axs[0,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp2[qx]), float(new_temp2[qx])), marker=None, linestyle="-", lw=3, c="orange", alpha=0.6)
        else:
            axs[0,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp1[qx]), float(new_temp1[qx])), marker=None, linestyle="-", lw=3, c="lightgray")
            axs[0,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp2[qx]), float(new_temp2[qx])), marker=None, linestyle="-", lw=3, c="lightgray", alpha=0.6)            
        qx=qx+1

    #Plot copy number seg for Lung Pl #3
    #axs[6,1].set_xlim([-1000000, 157040895])
    #axs[6,1].set_ylim([-0.1, 2.1]) 
    #axs[6,1].set_yticks([0, 1, 2]) 
    #axs[6,1].set_yticklabels([0, 1, 2])
    df_pur = df_abs_cn[df_abs_cn['contig']=="hs9"]
    df_pur = df_pur[df_pur['start']>=21700000]
    df_pur = df_pur[df_pur['start']<=21900000]
    df_pur = df_pur[['LungPl3_A','LungPl3_B']]
    vals = df_pur.values.tolist()
    sum_list = []
    for tempX in vals:
        sum_list.append(sum(tempX))
    normal_contam_estimate = sum(sum_list)/len(sum_list)
    purity = 1 - normal_contam_estimate

    df_af = df_abs_cn[df_abs_cn['contig']=="hsX"]
    df_af = df_af[df_af['start']>=2800000]
    df_af = df_af[df_af['start']<=155600000]
    temp1 = list(df_af['LungPl3_A'])
    temp2 = list(df_af['LungPl3_B'])

    ploidy = 1.99

    new_temp1 = []
    new_temp2 = []
    for qq in temp1:
        if qq >= 0:
            new_temp1.append(qq)
            #new_temp1.append((qq*(ploidy*purity+2*(1-purity))-(1-purity))/purity)
        else:
            new_temp1.append(0)
    
    for qq in temp2:
        if qq >= 0:
            new_temp2.append(qq)
            #new_temp2.append((qq*(ploidy*purity+2*(1-purity))-(1-purity))/purity)
        else:
            new_temp2.append(0)

    start_vals = df_af['start'].tolist()
    end_vals = df_af['end'].tolist()
    qx=0
    while qx<len(start_vals):
        #axs[6,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp1[qx]), float(new_temp1[qx])), marker=None, linestyle="-", lw=3, c="black")
        #axs[6,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp2[qx]), float(new_temp2[qx])), marker=None, linestyle="-", lw=3, c="orange", alpha=0.6)
        qx=qx+1

    df = input_counts_with_phase_df
    #Filter to sites with at least 10 RNA reads
    df = df[df['totalCount']>10]

    #Implement DNA read depth + germline VAF standard deviation filters
    df = df[df['position'].isin(valid_pos)]

    #Filter escapee exons + PARs + regions where allele A and allele B copy number are not 1
    valid_new_starts = valid_new_starts + invalid_starts
    valid_new_ends = valid_new_ends + invalid_ends

    invalid_pos = []
    all_positions = df['position'].tolist()
    w=0
    while w<len(valid_new_starts):
        for z in all_positions:
            if z >= int(valid_new_starts[w]) and z <= int(valid_new_ends[w]):
                invalid_pos.append(z)
        w=w+1

    df = df[~(df['position'].isin(invalid_pos))]
    df.to_csv("/Users/ananthansadagopan/Downloads/temp3.csv", index=False)

    #Calculate minor allele fraction
    df['AF'] = df['altCount']/(df['refCount']+df['altCount'])

    #Filter out points with really low DNA MAF in primary
    df['DNA_AF_primary'] = df['DNA_alt_count'].astype(int)/df['total_DNA_Count'].astype(int)
    df['DNA_Lung_total'] = df['DNA_ref_count_Lung'].astype(int) + df['DNA_alt_count_Lung'].astype(int)
    df['DNA_AF_lung'] = df['DNA_alt_count_Lung'].astype(int)/df['DNA_Lung_total'].astype(int)
    df = df[df['DNA_AF_primary']<=0.9]
    df = df[df['DNA_AF_primary']>=0.1]

    xs = df['position'].tolist()
    ys = df['AF'].tolist()
    MAF_vals = [min(x,1-x) for x in ys]
    df['MAF'] = MAF_vals
    phase_temp = df['Phase'].tolist()

    df.to_csv("/Users/ananthansadagopan/Downloads/livertum2.csv")
    
    """
    phase_temp = []
    hap_chr = hapbased_phase[hapbased_phase['chr']=="chrX"]
    hap_pos = hap_chr['pos'].tolist()
    final_hap = hap_chr['final_hap'].tolist()
    for t2pos in xs:
        wx=0
        while wx<len(hap_pos):
            if t2pos <= hap_pos[wx]:
                if final_hap[wx-1] == 1:
                    phase_temp.append("0|1")
                elif final_hap[wx-1] == -1:
                    phase_temp.append("1|0")
                else:
                    print("ERROR")
                break
            else:
                wx=wx+1
    """
    
    ys_DNA_primary = df['DNA_AF_primary'].tolist()
    ys_DNA_lung = df['DNA_AF_lung'].tolist()

    clist = []
    clist_inv = []

    a=0
    while a<len(phase_temp):
        if xs[a] <= 49028726: #49028726
            clist.append("lightgray")
            clist_inv.append("lightgray")
        else:
            if phase_temp[a] == "1|0":
                clist.append("black")
                clist_inv.append("orange")
            elif phase_temp[a] == "0|1":
                clist.append("orange")
                clist_inv.append("black")
        a=a+1
    #ys = [min(x,1-x) for x in ys_unp]
    #df['MAF'] = ys

    #Plot allele fraction vs. position

    axs[3,1].scatter(xs,ys, c=clist, alpha=1, s=6.5)
    axs[3,1].scatter(xs,[1-y for y in ys], c=clist_inv, alpha=1, s=6.5)

    axs[3,1].set_facecolor("white")
    axs[3,1].grid(False)
    axs[3,1].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4) 
    axs[3,1].spines['bottom'].set_color('0')
    axs[3,1].spines['left'].set_color('0')

    #Plot DNA allele fraction vs. position

    axs[2,1].scatter(xs,ys_DNA_primary, c=clist, alpha=1, s=6.5)
    axs[2,1].scatter(xs,[1-y for y in ys_DNA_primary], c=clist_inv, alpha=1, s=6.5)

    axs[2,1].set_facecolor("white")
    axs[2,1].grid(False)
    axs[2,1].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4) 
    axs[2,1].tick_params(axis='y', which='major', left=True, labelleft=True, labelsize=10, size=4) 
    axs[2,1].spines['bottom'].set_color('0')
    axs[2,1].spines['left'].set_color('0')

    #axs[6,1].scatter(xs,ys_DNA_lung, c=clist, alpha=1, s=6.5)
    #axs[6,1].scatter(xs,[1-y for y in ys_DNA_lung], c=clist_inv, alpha=1, s=6.5)

    #axs[6,1].set_facecolor("white")
    #axs[6,1].grid(False)
    #axs[6,1].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4) 
    #axs[6,1].tick_params(axis='y', which='major', left=True, labelleft=True, labelsize=10, size=4) 
    #axs[6,1].spines['bottom'].set_color('0')
    #axs[6,1].spines['left'].set_color('0')

    #Plot CDF of MAF (expressed SNVs + biallelic SNVs)
    count, bins_count = np.histogram(xs, bins=5000)
    pdf = count / sum(count)
    if reverse_val == True:
        cdf = 1-np.cumsum(pdf)
    else:
        cdf = np.cumsum(pdf)
    axs[4,1].plot(bins_count[1:], cdf, c="black")
    axs[4,1].set_xlim([-1000000, 157040895])
    axs[4,1].set_ylim([-0.01, 1.01])

    x_max = bins_count[np.argmax(cdf >= 1)]
    x_min = bins_count[np.argmax(cdf <= 0)]
    cdf_min = np.min(cdf)
    print(cdf_min)
    if reverse_val == True:
        axs[4,1].plot((0, x_max), (1, 1), linewidth=1.8, c="red", marker=None, dash_capstyle="butt")
        axs[4,1].plot((x_min, 157040895), (0, cdf_min), linewidth=1.8, c="red", marker=None, dash_capstyle="butt") 
    else:
        axs[4,1].plot((x_max, 157040895), (1, 1), linewidth=1.8, c="red", marker=None, dash_capstyle="butt")
        axs[4,1].plot((0, x_min-500000), (0, cdf_min), linewidth=1.8, c="red", marker=None, dash_capstyle="butt")
    sub_xs = []
    sub_ys = []

    xy=0
    while xy<len(xs):
        if ys[xy] >= 0.2 and ys[xy] <= 0.8:
            sub_xs.append(xs[xy])
            sub_ys.append(ys[xy])
        xy=xy+1

    count, bins_count = np.histogram(sub_xs, bins=5000)
    pdf = count / sum(count)

    if reverse_val == True:
        cdf = 1-np.cumsum(pdf)
    else:
        cdf = np.cumsum(pdf)
    axs[4,1].plot(bins_count[1:], cdf, c="red")
    axs[4,1].set_xlim([-1000000, 157040895])
    axs[4,1].set_ylim([-0.01, 1.01])

    x_max = bins_count[np.argmax(cdf >= 1)]
    x_min = bins_count[np.argmax(cdf <= 0)]
    cdf_min = np.min(cdf)
    print(cdf_min)
    if reverse_val == True:
        axs[4,1].plot((0, x_max), (1, 1), linewidth=1.8, c="red", marker=None, dash_capstyle="butt")
        axs[4,1].plot((x_min, 157040895), (0, cdf_min), linewidth=1.8, c="red", marker=None, dash_capstyle="butt") 
    else:
        axs[4,1].plot((x_max, 157040895), (1, 1), linewidth=1.8, c="red", marker=None, dash_capstyle="butt")
        axs[4,1].plot((0, x_min-500000), (0, cdf_min), linewidth=1.8, c="red", marker=None, dash_capstyle="butt")

    patch1 = mpatches.Patch(color="red", label="Biallelic\nSNVs (RNA)")
    patch2 = mpatches.Patch(color="black", label="Expressed\nSNVs")

    legend = axs[4,1].legend(handles=[patch1, patch2], bbox_to_anchor=(-0.25, 0.84))

    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')

    patch1 = mpatches.Patch(color="black", label="HapA (chrXa)")
    patch2 = mpatches.Patch(color="orange", label="HapB (chrXi)")

    legend = axs[3,1].legend(handles=[patch1, patch2], bbox_to_anchor=(-0.25, 0.64), title="Haplotype")

    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')

    #patch1 = mpatches.Patch(color="black", label="HapA (chrXa)")
    #patch2 = mpatches.Patch(color="orange", label="HapB (chrXi)")

    #legend = axs[2,1].legend(handles=[patch1, patch2], bbox_to_anchor=(-0.25, 1), title="Copy Number")

    #frame = legend.get_frame()
    #frame.set_facecolor('white')
    #frame.set_edgecolor('black')

    axs[4,1].set_facecolor("white")
    axs[4,1].grid(False)
    axs[4,1].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4) 
    axs[4,1].tick_params(axis='y', which='major', left=True, labelleft=True, labelsize=10, size=4) 
    axs[4,1].spines['bottom'].set_color('0')
    axs[4,1].spines['left'].set_color('0')

    #Plot MAF distributions left and right of break
    distal = df[(df['position']<49028726)]
    proximal = df[(df['position']>=49043410)]
    distal_vals = [x for x in distal['MAF'].tolist() if x==x]
    prox_vals = [x for x in proximal['MAF'].tolist() if x==x]

    median_maf_distal = statistics.median(distal_vals)
    median_maf_prox = statistics.median(prox_vals)

    try:
        wz,p = scipy.stats.ks_2samp(prox_vals, distal_vals, alternative="two-sided")
        wz, p = scipy.stats.mannwhitneyu(prox_vals,distal_vals, alternative="two-sided")
        KS_pval.append(p)
        axs[6,1].set_title('P = %s' % ("{:.2E}".format(p)), fontdict={'fontsize': 18}, y=-1.5)
    except:
        wz=0
        p=1
        KS_pval.append(p)

    axs[3,1].set_ylim([-0.01, 1.02])
    axs[2,1].set_ylim([-0.01, 1.02])
    #axs[6,1].set_ylim([-0.01, 1.02])

    #axs[5,1].axvspan(0, 49028726, alpha=0.15, color='red')
    #axs[5,1].axvspan(49043410, 156040895, alpha=0.15, color='blue')

    #Label axes
    axs[4,1].set_ylabel("Cumulative Probability")
    axs[4,1].set_xlabel("hg38 chrX Position (bp)", labelpad=5)
    axs[3,1].set_ylabel("Allele Fraction (RNA, LiverTum2)", labelpad=5)
    axs[2,1].set_ylabel("Allele Fraction (DNA, LiverTum2)", labelpad=5)
    #axs[6,1].set_ylabel("Allele Fraction (DNA, Lung Pl #1)", labelpad=5)
    axs[3,1].set_xlabel("hg38 chrX Position (bp)", labelpad=5)
    axs[2,1].set_xlabel("hg38 chrX Position (bp)", labelpad=5)
    #axs[1,2].set_xlabel("Density\nRight of TFE3")
    #axs[3,1].set_xlabel("hg38 chrX Position (bp)")
    #axs[2,1].set_xlabel("hg38 chrX Position (bp)") 
    #axs[6,1].set_ylabel("Copy\nNumber")
    #axs[6,1].set_xlabel("hg38 chrX Position (bp)")  
    #axs[5,1].yaxis.tick_right()
    axs[3,1].set_xlim([-1000000, 157040895])

    #Formatting
    axs[1,0].set_facecolor("white")
    axs[1,0].grid(False)

    axs[1,2].set_facecolor("white")
    axs[1,2].grid(False)


    axs[0,1].set_xlim([-1000000, 157040895])
    axs[1,1].set_xlim([-1000000, 157040895])

    axs[0,1].set_xlabel("hg38 chrX Position (bp)", labelpad=5)
    axs[1,1].set_xlabel("hg38 chrX Position (bp)", labelpad=5)
    axs[0,1].set_ylabel("Allelic\nDepth")
    axs[1,1].set_ylabel("Allelic\nDepth")
    axs[0,1].set_ylim([-0.1, 1.1])
    axs[1,1].set_ylim([-0.1, 1.1])
    axs[0,1].axvline(49028726, c="black", lw=0.5, ls="--")
    axs[1,1].axvline(49028726, c="black", lw=0.5, ls="--")

    axs[0,1].set_facecolor("white")
    axs[0,1].grid(False)
    axs[0,1].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4)
    axs[0,1].tick_params(axis='y', which='major', left=True, labelsize=10, size=4, labelleft=True)
    axs[0,1].spines['bottom'].set_color('0')
    axs[0,1].spines['left'].set_color('0')

    axs[1,1].set_facecolor("white")
    axs[1,1].grid(False)
    axs[1,1].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4)
    axs[1,1].tick_params(axis='y', which='major', left=True, labelsize=10, size=4, labelleft=True)
    axs[1,1].spines['bottom'].set_color('0')
    axs[1,1].spines['left'].set_color('0')

    axs[5,2].get_xaxis().set_visible(False)
    axs[5,2].get_yaxis().set_visible(False)
    axs[5,0].get_xaxis().set_visible(False)
    axs[5,0].get_yaxis().set_visible(False)

    axs[5,2].set_facecolor("white")
    axs[5,2].grid(False)

    axs[5,0].set_facecolor("white")
    axs[5,0].grid(False)

    axs[6,2].get_xaxis().set_visible(False)
    axs[6,2].get_yaxis().set_visible(False)
    axs[6,0].get_xaxis().set_visible(False)
    axs[6,0].get_yaxis().set_visible(False)

    axs[6,2].set_facecolor("white")
    axs[6,2].grid(False)

    axs[6,0].set_facecolor("white")
    axs[6,0].grid(False)


    """
    axs[3,1].set_facecolor("white")
    axs[3,1].grid(False)
    axs[3,1].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4)
    axs[3,1].tick_params(axis='y', which='major', left=True, labelsize=10, size=4, labelleft=True)
    axs[3,1].spines['bottom'].set_color('0')
    axs[3,1].spines['left'].set_color('0')

    axs[2,1].set_facecolor("white")
    axs[2,1].grid(False)
    axs[2,1].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4)
    axs[2,1].tick_params(axis='y', which='major', left=True, labelsize=10, size=4, labelleft=True)
    axs[2,1].spines['bottom'].set_color('0')
    axs[2,1].spines['left'].set_color('0')

    axs[6,1].set_facecolor("white")
    axs[6,1].grid(False)
    axs[6,1].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4)
    axs[6,1].tick_params(axis='y', which='major', left=True, labelsize=10, size=4, labelleft=True)
    axs[6,1].spines['bottom'].set_color('0')
    axs[6,1].spines['left'].set_color('0')
    """

    #axs[4,1].spines['bottom'].set_color('0')
    #axs[4,1].spines['left'].set_color('0')

    abc=2
    while abc<=4:

        axs[abc,2].set_facecolor("white")
        axs[abc,2].grid(False)

        axs[abc,0].set_facecolor("white")
        axs[abc,0].grid(False)

        axs[abc,2].get_xaxis().set_visible(False)
        axs[abc,2].get_yaxis().set_visible(False)
        axs[abc,0].get_xaxis().set_visible(False)
        axs[abc,0].get_yaxis().set_visible(False)

        abc=abc+1

    axs[0,2].set_facecolor("white")
    axs[0,2].grid(False)

    axs[0,0].set_facecolor("white")
    axs[0,0].grid(False)

    axs[0,2].get_xaxis().set_visible(False)
    axs[0,2].get_yaxis().set_visible(False)
    axs[0,0].get_xaxis().set_visible(False)
    axs[0,0].get_yaxis().set_visible(False)

    """
    axs[3,2].set_facecolor("white")
    axs[3,2].grid(False)

    axs[3,0].set_facecolor("white")
    axs[3,0].grid(False)

    axs[3,2].get_xaxis().set_visible(False)
    axs[3,2].get_yaxis().set_visible(False)
    axs[3,0].get_xaxis().set_visible(False)
    axs[3,0].get_yaxis().set_visible(False)

    axs[4,2].set_facecolor("white")
    axs[4,2].grid(False)

    axs[4,0].set_facecolor("white")
    axs[4,0].grid(False)

    axs[4,2].get_xaxis().set_visible(False)
    axs[4,2].get_yaxis().set_visible(False)
    axs[4,0].get_xaxis().set_visible(False)
    axs[4,0].get_yaxis().set_visible(False)

    axs[5,2].set_facecolor("white")
    axs[5,2].grid(False)

    axs[5,0].set_facecolor("white")
    axs[5,0].grid(False)

    axs[5,2].get_xaxis().set_visible(False)
    axs[5,2].get_yaxis().set_visible(False)
    axs[5,0].get_xaxis().set_visible(False)
    axs[5,0].get_yaxis().set_visible(False)
    """

    axs[1,2].get_xaxis().set_visible(False)
    axs[1,2].get_yaxis().set_visible(False)
    axs[1,0].get_xaxis().set_visible(False)
    axs[1,0].get_yaxis().set_visible(False)

    axs[4,1].axvline(49028726, c="black", lw=0.5, ls="--")
    axs[3,1].axvline(49028726, c="black", lw=0.5, ls="--")
    axs[2,1].axvline(49028726, c="black", lw=0.5, ls="--")
    #axs[3,1].axvline(49028726, c="black", lw=0.5, ls="--")
    #axs[2,1].axvline(49028726, c="black", lw=0.5, ls="--")
    #axs[6,1].axvline(49028726, c="black", lw=0.5, ls="--")

    #axs[5,1].axhline(median_maf_distal, 1000000/158040895, 50028726/158040895, c="blue", lw=0.5, ls="--")
    #axs[5,1].axhline(median_maf_prox, 50043410/158040895, 157040895/158040895, c="blue", lw=0.5, ls="--")

    axs[3,1].tick_params(axis='y', which='major', left=True, labelsize=10, size=4, labelleft=True)

    fig.tight_layout()
    plt.subplots_adjust(wspace=0.55)
    plt.subplots_adjust(hspace=0.44)

    fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/Xi_Xa_fusion/BTRCC18_phased/" + www + "_AF_vs_position_chrX_no_CN_LiverTum2.pdf", dpi=dpi_set, bbox_inches = 'tight')

df_out = pd.DataFrame([vals, KS_pval]).T
df_out.columns = ['Sample', 'pval_KS']
df_out.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/Xi_Xa_fusion/BTRCC18_phased/BTRCC18_XR_summary_KS_test_NEW_paired_LiverTum2.csv")
