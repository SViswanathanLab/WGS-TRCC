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
sns.set(rc={'figure.figsize':(5.5,6)})
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
df_abs_cn = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/WGS_merged_seg_file_placeholder.csv")

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

#Specify directory with ASEreadcounter out files
temp_path = "/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/RNAseq_tophat_paired/"
vals = [f for f in listdir(temp_path) if "chrX_RNA_paired_on_germline_hets_unphased" in f and "13" not in f and "18" in f]
KS_pval = []

#Iterate over each file
for www in vals:
    value = www.split("chrX_RNA_paired_on_germline_hets_unphased_")[1].split(".log")[0]
    temp_chr = "chrX"
    temp_id = str(value.split("CC_")[-1].split("CC")[-1].split("_")[0])
    if temp_id == str(8):
        temp_id = "DTRCC8"

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
        print(df.columns.tolist())
        if "17" in a:
            temp_val = df["RTRCC_17_Blood_Normal"].tolist()
        else:
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

    temp_file = temp_path + www
    df = pd.read_csv(temp_file, sep="\t")

    #Filter to sites with at least 10 RNA reads
    df = df[df['totalCount']>10]

    #Implement DNA read depth + germline VAF standard deviation filters
    df = df[df['position'].isin(valid_pos)]





    if "18" in str(temp_id):
        print("TRCC18_CONFIRMED")

        xs = df['position'].tolist()

        df_abs_cn = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/BTRCC18_allelic_CN.xlsx")

        df_pur = df_abs_cn[df_abs_cn['contig']=="hs9"]
        df_pur = df_pur[df_pur['start']>=21700000]
        df_pur = df_pur[df_pur['start']<=21900000]
        df_pur = df_pur[['Primary_A','Primary_B']]
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
        temp1 = list(df_af['Primary_A'])
        temp2 = list(df_af['Primary_B'])

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
                if start_vals[qx]>49028726:
                    pass
                else:
                    pass
            else:
                pass
            if new_temp1[qx] < meantemp1-stdtemp1 or new_temp1[qx] > meantemp1+stdtemp1 or new_temp2[qx] < meantemp2-stdtemp2 or new_temp2[qx] > meantemp2+stdtemp2: #Identify sites with CN imbalance
                invalid_starts.append(start_vals[qx])
                invalid_ends.append(end_vals[qx])
            qx=qx+1
            valid_new_starts = valid_new_starts + invalid_starts
            valid_new_ends = valid_new_ends + invalid_ends



    #Filter escapee exons + PARs
    invalid_pos = []
    all_positions = df['position'].tolist()
    w=0
    while w<len(valid_new_starts):
        for z in all_positions:
            if z >= int(valid_new_starts[w]) and z <= int(valid_new_ends[w]):
                invalid_pos.append(z)
        w=w+1

    df = df[~(df['position'].isin(invalid_pos))]

    #Calculate minor allele fraction
    df['AF'] = df['refCount']/(df['refCount']+df['altCount'])
    xs = df['position'].tolist()
    ys_unp = df['AF'].tolist()
    ys = [min(x,1-x) for x in ys_unp]
    df['MAF'] = ys

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

    fig, axs = plt.subplots(2, 1, sharey='row', gridspec_kw={'height_ratios': [6, 1]})
    d = {}
    yticks = []
    yticklabels = []

    # ideogram.txt downloaded from UCSC's table browser
    for xranges, yrange, colors, label in ideograms('/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/cytoBandIdeo_chrX.txt'):
        print(xranges)
        coll = BrokenBarHCollection(xranges, yrange, facecolors=colors, lw=0)
        axs[1].add_collection(coll)
        center = yrange[0] + yrange[1]/2.
        yticks.append(center)
        yticklabels.append(label)
        d[label] = xranges

    axs[1].axis('tight')
    axs[1].set_yticks([])
    axs[1].set_yticklabels([])
    axs[1].set_xticks([])

    axs[1].set_xlim([-1000000, 157040895])

    #Plot CDF of MAF (expressed SNVs + biallelic SNVs)
    count, bins_count = np.histogram(xs, bins=5000)
    pdf = count / sum(count)
    cdf = np.cumsum(pdf)
    axs[0].plot(bins_count[1:], cdf, c="black")
    axs[0].set_xlim([-1000000, 157040895])
    axs[0].set_ylim([-0.01, 1.01])

    x_max = bins_count[np.argmax(cdf >= 1)]
    x_min = bins_count[np.argmax(cdf <= 0)]
    cdf_min = np.min(cdf)
    print(x_max)
    axs[0].plot((max(xs), 157040895), (1, 1), linewidth=1.4, c="black", marker=None, dash_capstyle="butt")
    axs[0].plot((0, min(xs)-500000), (0, cdf_min), linewidth=1.4, c="black", marker=None, dash_capstyle="butt")

    sub_xs = []
    sub_ys = []

    xy=0
    while xy<len(xs):
        if ys[xy] >= 0.2:
            sub_xs.append(xs[xy])
            sub_ys.append(ys[xy])
        xy=xy+1

    count, bins_count = np.histogram(sub_xs, bins=5000)
    pdf = count / sum(count)
    cdf = np.cumsum(pdf)
    axs[0].plot(bins_count[1:], cdf, c="red")
    axs[0].set_xlim([-1000000, 157040895])
    axs[0].set_ylim([-0.01, 1.01])

    x_max = bins_count[np.argmax(cdf >= 1)]
    x_min = bins_count[np.argmax(cdf <= 0)]
    cdf_min = np.min(cdf)
    print(cdf_min)
    axs[0].plot((max(sub_xs), 157040895), (1, 1), linewidth=1.4, c="red", marker=None, dash_capstyle="butt")
    axs[0].plot((0, min(sub_xs)-500000), (0, cdf_min), linewidth=1.4, c="red", marker=None, dash_capstyle="butt")

    #patch1 = mpatches.Patch(color="red", label="Biallelic\nSNVs (RNA)")
    #patch2 = mpatches.Patch(color="black", label="Expressed\nSNVs")

    #legend = axs[0,1].legend(handles=[patch1, patch2], bbox_to_anchor=(-0.25, 0.84))

    #frame = legend.get_frame()
    #frame.set_facecolor('white')
    #frame.set_edgecolor('black')

    axs[0].set_facecolor("white")
    axs[0].grid(False)
    axs[0].tick_params(axis='x', which='major', bottom=True, labelsize=16, size=4) 
    axs[0].tick_params(axis='y', which='major', left=True, labelleft=True, labelsize=16, size=4) 
    axs[0].spines['bottom'].set_color('0')
    axs[0].spines['left'].set_color('0')

    axs[0].set_xlabel("hg38 chrX Position (bp)", labelpad=5, size=20)
    axs[0].set_ylabel("Cumulative Probability", labelpad=5, size=20)
    axs[0].set_xlim([-1000000, 157040895])

    axs[1].set_facecolor("white")
    axs[1].grid(False)

    axs[0].axvline(49028726, c="black", lw=1.5, ls="--")

    axs[1].tick_params(axis='y', which='major', left=True, right=False, labelsize=16, size=4, labelleft=True)

    fig.tight_layout()
    plt.subplots_adjust(wspace=1)

    fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/RNAseq_tophat_paired_plots/TRUNCATED_" + www + "_AF_vs_position.pdf", dpi=dpi_set, bbox_inches = 'tight')

#df_out = pd.DataFrame([vals, KS_pval]).T
#df_out.columns = ['Sample', 'pval_KS']
#df_out.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/WGS/WGS_XR_summary_KS_test_NEW_paired.csv")