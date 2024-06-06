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

#Compiled seg file
df_abs_cn = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/BTRCC18_allelic_CN.xlsx")

valid_new_starts = []
valid_new_ends = []
KS_pval = []
#lengths = [135086622, 248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468]
#main_chr_list = ["chr11", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
#lengths = [133797422]
#main_chr_list = ["chr10"]
lengths = [159345973]
main_chr_list = ["chr7"]

vals = []

for a in lengths:
    vals.append("BTRCC18")

qqx=0

print("reading VCF")
all_files_in_dir = os.listdir("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/germline_hets/")
temp_DNA_id = [x for x in all_files_in_dir if '18' in x]
path1 = "/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/germline_hets/" + temp_DNA_id[0]
#names, header_temp = get_vcf_names(path1)
#dna_df = pd.read_csv(path1, comment='#', delim_whitespace=True, header=None, names=names, compression="gzip")

path2 = "/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/Xi_Xa_fusion/BTRCC18_phased/chr7_het_eagle2phased.vcf.gz"
names, header_temp = get_vcf_names(path2)
df = pd.read_csv(path2, comment='#', delim_whitespace=True, header=None, names=names, compression="gzip")
print("read VCF")

tumors_with_Xq_del = ['LeftMedLNTum_1','LiverTum_4']
normals = ['B_TRCC_18_Normal', 'KidneyNormal', 'LiverNormal']

ref_df_rows = []
alt_df_rows = []

for a in tumors_with_Xq_del:
    temp_val = df[a].tolist()
    new_val1 = []
    new_val2 = []
    for q in temp_val:
        temp = q.split(":")[1]
        temp2 = q.split(":")[0]
        sub = temp.split(",")
        if sub != ['.']:
            new_val1.append(sub[0])
            new_val2.append(sub[1])
        else:
            new_val1.append(float("nan"))
            new_val2.append(float("nan"))          
    ref_df_rows.append(new_val1)
    alt_df_rows.append(new_val2)

ref_df = pd.DataFrame(ref_df_rows).T
alt_df = pd.DataFrame(alt_df_rows).T
ref_df_rows = ref_df.values.tolist()
alt_df_rows = alt_df.values.tolist()

AF = []

q=0
while q<len(ref_df_rows):
    temp_row1 = [int(x) for x in ref_df_rows[q] if x == x]
    temp_row2 = [int(x) for x in alt_df_rows[q] if x == x]
    try:
        sum1 = sum(temp_row1)/len(temp_row1)
        sum2 = sum(temp_row2)/len(temp_row2)
        AF.append(sum1/(sum1 + sum2))
    except:
        AF.append(float("nan"))
    q=q+1
    
df['AF_Xq_del'] = AF

ref_df_rows = []
alt_df_rows = []

for a in normals:
    temp_val = df[a].tolist()
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

ref_df = pd.DataFrame(ref_df_rows).T
alt_df = pd.DataFrame(alt_df_rows).T
ref_df_rows = ref_df.values.tolist()
alt_df_rows = alt_df.values.tolist()

AF = []

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
    else:
        AF.append(sum1/(sum1 + sum2))
    q=q+1
    
df['AF_normal'] = AF
df = df[df['AF_normal'].notna()]
df = df[df['AF_Xq_del'].notna()]

AF_n_vals = df['AF_normal'].tolist()
AF_Xq_vals = df['AF_Xq_del'].tolist()

phase = []

a=0
while a<len(AF_n_vals):
    if AF_n_vals[a] < AF_Xq_vals[a]:
        phase.append("0|1")
    else:
        phase.append("1|0")
    a=a+1

df['Phase'] = phase
phase_dict = dict(zip(df['POS'].tolist(), df['Phase'].tolist()))

dna_df = df
















#Iterate over each file
for www in vals:

    temp_chr = main_chr_list[qqx]
    hg38_chr_length = lengths[qqx]

    value = "B_TRCC_18_Tumor"
    temp_id = '18'
    hs_val = "hs" + str(temp_chr.split("chr")[-1])

    #Find germline hets with DNA VAF within 1SD of the mean across chr
    sd_cutoff = 1
    df = dna_df[dna_df['#CHROM']==temp_chr]
    
    normals = [temp_id]
    ref_df_rows = []
    alt_df_rows = []

    for a in normals:
        temp_val = df["B_TRCC_18_Normal"].tolist()
        new_val1 = []
        new_val2 = []
        for q in temp_val:
            temp = q.split(":")[1]
            temp2 = q.split(":")[0]
            sub = temp.split(",")
            if sub != ['.'] and temp2 != "0/0" and temp2 != "1/1" and temp2 != "0|0" and temp2 != "1|1":
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
    
    #Repeat for Tumor
    normals = ['B_TRCC_18_Tumor']
    ref_df_rows = []
    alt_df_rows = []

    for a in normals:
        temp_val = df[a].tolist()
        new_val1 = []
        new_val2 = []
        for q in temp_val:
            temp = q.split(":")[1]
            temp2 = q.split(":")[0]
            sub = temp.split(",")
            if sub != ['.']:# and temp2 != "0/0" and temp2 != "1/1":
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
        
    df['total_tumor'] = total_sum
    df['AF_tumor'] = AF

    df = df[df['AF_normal'].notna()]
    df = df[df['total_sum'].notna()]

    tumor_AF_vals = df['AF_tumor'].tolist()
    tumor_pos_vals = df['POS'].tolist()

    pos_tumor_AF_dict = dict(zip(tumor_pos_vals,tumor_AF_vals))

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
        xranges, colors = [(0,3000000)], [(.7, .7, .7)]
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
    d = {}
    yticks = []
    yticklabels = []

    temp_ideo = '/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/cytoBandIdeo.txt'
    ideo_df = pd.read_csv(temp_ideo, sep="\t", header=None)
    ideo_df = ideo_df[ideo_df[0]==temp_chr]
    ideo_df.to_csv('/Users/ananthansadagopan/Downloads/temp_ideo.txt', sep="\t", index=False, header=False)

    # ideogram.txt downloaded from UCSC's table browser
    for xranges, yrange, colors, label in ideograms('/Users/ananthansadagopan/Downloads/temp_ideo.txt'):
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

    axs[5,1].set_xlim([-1000000, hg38_chr_length+1000000])

    # ideogram.txt downloaded from UCSC's table browser
    for xranges, yrange, colors, label in ideograms('/Users/ananthansadagopan/Downloads/temp_ideo.txt'):
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

    axs[6,1].set_xlim([-1000000, hg38_chr_length+1000000])

    #Plot copy number seg for Primary
    #axs[3,1].set_xlim([-1000000, hg38_chr_length+1000000])
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

    df_af = df_abs_cn[df_abs_cn['contig']==hs_val]
    temp1 = list(df_af['LiverTum2_A'])
    temp2 = list(df_af['LiverTum2_B'])

    ploidy = 2.04

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
            axs[1,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp1[qx]), float(new_temp1[qx])), marker=None, linestyle="-", lw=3, c="black")
            axs[1,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp2[qx]), float(new_temp2[qx])), marker=None, linestyle="-", lw=3, c="orange", alpha=0.6)
        else:
            #pass
            axs[1,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp1[qx]), float(new_temp1[qx])), marker=None, linestyle="-", lw=3, c="green", zorder=10)
            axs[1,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp2[qx]), float(new_temp2[qx])), marker=None, linestyle="-", lw=3, c="green", zorder=10)  
        if new_temp1[qx] < meantemp1-stdtemp1 or new_temp1[qx] > meantemp1+stdtemp1 or new_temp2[qx] < meantemp2-stdtemp2 or new_temp2[qx] > meantemp2+stdtemp2: #Identify sites with CN imbalance
            invalid_starts.append(start_vals[qx])
            invalid_ends.append(end_vals[qx])
        qx=qx+1

    #Plot copy number seg for Lung Pl #1
    #axs[2,1].set_xlim([-1000000, hg38_chr_length+1000000])
    #axs[2,1].set_ylim([-0.1, 2.1]) 
    #axs[2,1].set_yticks([0, 1, 2]) 
    #axs[2,1].set_yticklabels([0, 1, 2]) 

    df_pur = df_abs_cn[df_abs_cn['contig']=="hs9"]
    df_pur = df_pur[df_pur['start']>=21700000]
    df_pur = df_pur[df_pur['start']<=21900000]
    df_pur = df_pur[['LeftMedLN1_A','LeftMedLN1_B']]
    vals = df_pur.values.tolist()
    sum_list = []
    for tempX in vals:
        sum_list.append(sum(tempX))
    normal_contam_estimate = sum(sum_list)/len(sum_list)
    #purity = 1 - normal_contam_estimate
    purity = 0.5475

    df_af = df_abs_cn[df_abs_cn['contig']==hs_val]
    temp1 = list(df_af['LeftMedLN1_A'])
    temp2 = list(df_af['LeftMedLN1_B'])

    ploidy = 3.747

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
        pass
        axs[0,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp1[qx]), float(new_temp1[qx])), marker=None, linestyle="-", lw=3, c="orange")
        axs[0,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp2[qx]), float(new_temp2[qx])), marker=None, linestyle="-", lw=3, c="black", alpha=0.6)
        qx=qx+1

    #Plot copy number seg for Lung Pl #3
    #axs[6,1].set_xlim([-1000000, hg38_chr_length+1000000])
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

    df_af = df_abs_cn[df_abs_cn['contig']==hs_val]
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
        pass
        #axs[6,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp1[qx]), float(new_temp1[qx])), marker=None, linestyle="-", lw=3, c="black")
        #axs[6,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(new_temp2[qx]), float(new_temp2[qx])), marker=None, linestyle="-", lw=3, c="orange", alpha=0.6)
        qx=qx+1

    ase_output = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/RNAseq_paired_all_chr/allchr_RNA_paired_on_germline_hets_unphased_LiverTum_2.log", sep="\t")
    df = ase_output[ase_output['contig']==temp_chr]
    #Filter to sites with at least 10 RNA reads
    df = df[df['totalCount']>10]

    #Implement DNA read depth + germline VAF standard deviation filters
    df = df[df['position'].isin(valid_pos)]

    df.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/Xi_Xa_fusion/BTRCC18_phased/chr7_unphased_DNA_and_RNA_LiverTum2.csv", index=False)

    """
    #Filter regions where allele A and allele B copy number are not 1
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
    """

    #Calculate minor allele fraction
    df['AF'] = df['altCount']/(df['refCount']+df['altCount'])
    print(df)

    xs = df['position'].tolist()
    ys = df['AF'].tolist()
    MAF_vals = [min(x,1-x) for x in ys]
    df['MAF'] = MAF_vals
    clist = []
    clist_inv = []

    a=0
    while a<len(xs):
        if xs[a] > hg38_chr_length: #49028726
            clist.append("lightgray")
            clist_inv.append("lightgray")
        else:
            if phase_dict[xs[a]] == "1|0":
                clist.append("black")
                clist_inv.append("orange")
            elif phase_dict[xs[a]] == "0|1":
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
    axs[3,1].tick_params(axis='y', which='major', left=True, labelleft=True, labelsize=10, size=4) 
    axs[3,1].spines['bottom'].set_color('0')
    axs[3,1].spines['left'].set_color('0')

    new_ys_vals = []
    for abcdefg in xs:
        try:
            new_ys_vals.append(float(pos_tumor_AF_dict[int(abcdefg)]))
        except:
            new_ys_vals.append(float("nan"))

    axs[2,1].scatter(xs,new_ys_vals, c=clist, alpha=1, s=6.5)
    axs[2,1].scatter(xs,[1-y for y in new_ys_vals], c=clist_inv, alpha=1, s=6.5)

    axs[2,1].set_facecolor("white")
    axs[2,1].grid(False)
    axs[2,1].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4) 
    axs[2,1].tick_params(axis='y', which='major', left=True, labelleft=True, labelsize=10, size=4) 
    axs[2,1].spines['bottom'].set_color('0')
    axs[2,1].spines['left'].set_color('0')

    #Plot CDF of MAF (expressed SNVs + biallelic SNVs)
    count, bins_count = np.histogram(xs, bins=5000)
    pdf = count / sum(count)
    cdf = np.cumsum(pdf)
    axs[4,1].plot(bins_count[1:], cdf, c="black")
    axs[4,1].set_xlim([-1000000, hg38_chr_length+1000000])
    axs[4,1].set_ylim([-0.01, 1])

    x_max = bins_count[np.argmax(cdf >= 0.9999)]
    x_min = bins_count[np.argmax(cdf <= 0)]
    cdf_min = np.min(cdf)
    print(cdf_min)
    print(max(xs))

    axs[4,1].plot((max(xs), hg38_chr_length), (1, 1), linewidth=1.4, c="black", marker=None, dash_capstyle="butt")
    axs[4,1].plot((0, x_min-500000), (0, cdf_min), linewidth=1.4, c="black", marker=None, dash_capstyle="butt")

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
    cdf = np.cumsum(pdf)
    axs[4,1].plot(bins_count[1:], cdf, c="red")
    axs[4,1].set_xlim([-1000000, hg38_chr_length+1000000])
    axs[4,1].set_ylim([-0.01, 1.02])

    x_max = bins_count[np.argmax(cdf >= 1)]
    x_min = bins_count[np.argmax(cdf <= 0)]
    cdf_min = np.min(cdf)
    axs[4,1].plot((x_max, hg38_chr_length), (1, 1), linewidth=1.4, c="red", marker=None, dash_capstyle="butt")
    axs[4,1].plot((0, x_min-500000), (0, cdf_min), linewidth=1.4, c="red", marker=None, dash_capstyle="butt")

    patch1 = mpatches.Patch(color="red", label="Biallelic\nSNVs (RNA)")
    patch2 = mpatches.Patch(color="black", label="Expressed\nSNVs")

    legend = axs[4,1].legend(handles=[patch1, patch2], bbox_to_anchor=(-0.25, 0.84))

    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')

    #patch1 = mpatches.Patch(color="black", label="HapA (chrXa)")
    #patch2 = mpatches.Patch(color="orange", label="HapB (chrXi)")

    #legend = axs[5,1].legend(handles=[patch1, patch2], bbox_to_anchor=(-0.25, 0.64), title="Haplotype")

    #frame = legend.get_frame()
    #frame.set_facecolor('white')
    #frame.set_edgecolor('black')

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
    distal = df[(df['position']<35176378)]
    proximal = df[(df['position']>=35176378)]
    distal_vals = [x for x in distal['MAF'].tolist() if x==x]
    prox_vals = [x for x in proximal['MAF'].tolist() if x==x]

    maf_vals_finalized = [float(x) for x in df['MAF'].tolist() if x == x]

    try:
        wz,p = scipy.stats.ks_2samp(prox_vals, distal_vals, alternative="two-sided")
        KS_pval.append(p)
        axs[6,1].set_title('P = %s' % ("{:.2E}".format(p)), fontdict={'fontsize': 18}, y=-1.5)
    except:
        wz=0
        p=1
        KS_pval.append(p)

    #axs[5,1].set_ylim([-0.01, 1.02])
    #axs[5,1].axvspan(0, 49028726, alpha=0.15, color='red')
    #axs[5,1].axvspan(49043410, 156040895, alpha=0.15, color='blue')

    #Label axes
    axs[4,1].set_ylabel("Cumulative Probability")
    axs[4,1].set_xlabel("hg38 " + temp_chr + " Position (bp)", labelpad=5)
    axs[3,1].set_ylabel("Allele Fraction (RNA)", labelpad=5)
    axs[3,1].set_xlabel("hg38 " + temp_chr + " Position (bp)", labelpad=5)
    axs[2,1].set_ylabel("Allele Fraction (DNA)", labelpad=5)
    axs[2,1].set_xlabel("hg38 " + temp_chr + " Position (bp)", labelpad=5)
    #axs[1,0].set_xlabel("Density\nLeft of TFE3")
    #axs[1,2].set_xlabel("Density\nRight of TFE3")
    #axs[3,1].set_ylabel("Copy\nNumber")
    #axs[3,1].set_xlabel("hg38 " + temp_chr + " Position (bp)")
    #axs[2,1].set_ylabel("Copy\nNumber")
    #axs[2,1].set_xlabel("hg38 " + temp_chr + " Position (bp)") 
    #xs[4,1].set_ylabel("Copy\nNumber")
    #axs[6,1].set_xlabel("hg38 " + temp_chr + " Position (bp)")  
    #axs[5,1].yaxis.tick_right()
    axs[4,1].set_xlim([-1000000, hg38_chr_length+1000000])
    axs[5,1].set_xlim([-1000000, hg38_chr_length+1000000])
    axs[3,1].set_xlim([-1000000, hg38_chr_length+1000000])
    axs[2,1].set_xlim([-1000000, hg38_chr_length+1000000])

    #Formatting
    axs[1,0].set_facecolor("white")
    axs[1,0].grid(False)

    axs[1,2].set_facecolor("white")
    axs[1,2].grid(False)


    axs[0,1].set_xlim([-1000000, hg38_chr_length+1000000])
    axs[1,1].set_xlim([-1000000, hg38_chr_length+1000000])

    axs[0,1].set_xlabel("hg38 " + temp_chr + " Position (bp)", labelpad=5)
    axs[1,1].set_xlabel("hg38 " + temp_chr + " Position (bp)", labelpad=5)
    axs[0,1].set_ylabel("Allelic\nDepth")
    axs[1,1].set_ylabel("Allelic\nDepth")
    axs[0,1].set_ylim([-0.1, 1.1])
    axs[1,1].set_ylim([-0.1, 1.1])

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

    axs[2,2].set_facecolor("white")
    axs[2,2].grid(False)

    axs[2,0].set_facecolor("white")
    axs[2,0].grid(False)

    axs[2,2].get_xaxis().set_visible(False)
    axs[2,2].get_yaxis().set_visible(False)
    axs[2,0].get_xaxis().set_visible(False)
    axs[2,0].get_yaxis().set_visible(False)

    axs[0,2].set_facecolor("white")
    axs[0,2].grid(False)

    axs[0,0].set_facecolor("white")
    axs[0,0].grid(False)

    axs[0,2].get_xaxis().set_visible(False)
    axs[0,2].get_yaxis().set_visible(False)
    axs[0,0].get_xaxis().set_visible(False)
    axs[0,0].get_yaxis().set_visible(False)

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
    axs[3,1].set_ylim([-0.01, 1.02])
    axs[2,1].set_ylim([-0.01, 1.02])

    #axs[4,1].axvline(35176378, c="black", lw=0.5, ls="--")
    #axs[5,1].axvline(49028726, c="black", lw=0.5, ls="--")
    #axs[3,1].axvline(35176378, c="black", lw=0.5, ls="--")
    #axs[2,1].axvline(35176378, c="black", lw=0.5, ls="--")
    #axs[6,1].axvline(49028726, c="black", lw=0.5, ls="--")

    fig.tight_layout()
    plt.subplots_adjust(wspace=0.55)
    plt.subplots_adjust(hspace=0.44)

    axs[5,1].tick_params(axis='y', which='major', left=True, labelsize=10, size=4, labelleft=True)
    #axs[5,1].axhline(statistics.median(maf_vals_finalized), c="blue", lw=0.5, ls="--")

    fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/Xi_Xa_fusion/BTRCC18_phased/" + www + "_AF_vs_position_" + temp_chr + "_no_CN_LiverTum2.pdf", dpi=dpi_set, bbox_inches = 'tight')
    qqx=qqx+1

#df_out = pd.DataFrame([vals, KS_pval]).T
#df_out.columns = ['Sample', 'pval_KS']
#df_out.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/Xi_Xa_fusion/BTRCC18_phased/BTRCC18_XR_summary_KS_test_NEW_paired_" + temp_chr + ".csv")
