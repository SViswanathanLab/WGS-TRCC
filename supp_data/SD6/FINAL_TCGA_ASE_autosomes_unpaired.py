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
sns.set(rc={'figure.figsize':(6.5,8.5)})
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

#TCGA ABSOLUTE CN
df_abs_cn = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/titan_calls/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg", sep="\t")









"""
#Defining escapee exon boundaries (creating list of excluded regions)
chrX_exons_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/chrX_all_exons_gencode_v43_UCSC.csv")
gene_classes_df = pd.read_excel("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/kidney_specific_chrX_gene_class_list_cotton.xlsx")
escape_gene = [x for x in gene_classes_df['Escape_plus_Variable'].tolist() if x == x]
#escape_gene = escape_gene + ['ACE2','ARSL','PUDP','TBL1X','TMSB4X','USP11']
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
"""

valid_new_starts = []
valid_new_ends = []

#Specify directory with ASEreadcounter out files
temp_path = "/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/TCGA/unpaired_all_chr_all_samples/"
vals = [f for f in listdir(temp_path) if "all_chr_TopHat2_RNA_unpaired_on_germline_hets_unphased_" in f]
KS_pval = []
distal_vals_total = []
prox_vals_total = []

id_list = ['TCGA-J7-8537', 'TCGA-B0-5705','TCGA-B8-5546','TCGA-BQ-7050','TCGA-G7-7501','TCGA-CJ-5681','TCGA-BP-4756','TCGA-6D-AA2E','TCGA-2Z-A9JO']
fusion_partner = ['DVL2','SFPQ','SFPQ','PRCC','SFPQ','KHSRP','SFPQ','MED15','U2AF2']
fusion_chromosome = ["chr17","chr1","chr1","chr1","chr1","chr19","chr1","chr22","chr17"]
fusion_location = [7128660,35176378,35176378,156800815,35176378,6413359,35176378,20587632,55674716]
temp_length = [83257441,248956422,248956422,248956422,248956422,58617616,248956422,50818468,83257441]

fusion_partner_dict = dict(zip(id_list, fusion_partner))
fusion_chromosome_dict = dict(zip(id_list, fusion_chromosome))
fusion_location_dict = dict(zip(id_list, fusion_location))
fusion_length_dict = dict(zip(id_list, temp_length))

#Iterate over each file
for www in vals:
    value = www.split("all_chr_TopHat2_RNA_unpaired_on_germline_hets_unphased_")[1].split(".log")[0]
    try:
        temp_id = value.split("_")[0]
    except:
        temp_id = value

    temp_chr=fusion_chromosome_dict[temp_id]
    temp_loc=fusion_location_dict[temp_id]
    temp_partner=fusion_partner_dict[temp_id]
    temp_length=fusion_length_dict[temp_id]

    print(temp_id)
    print(temp_chr)

    #Find germline hets with DNA VAF within 1SD of the mean across chrX
    sd_cutoff = 1
    print("reading VCF")
    path1 = "/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/TCGA/TCGA_VCFs/" + temp_id + "_germline_DNA.vcf.gz"
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

    temp_file = temp_path + temp_chr + www[7:]
    df = pd.read_csv(temp_file, sep="\t")
    df = df[df['contig']==temp_chr]

    print(df)

    #Filter to sites with at least 10 RNA reads
    df = df[df['totalCount']>10]

    print("TotalCount_Filtered")
    print(df)

    #Implement DNA read depth + germline VAF standard deviation filters
    df = df[df['position'].isin(valid_pos)]

    print("Position_Filtered")
    print(df)

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

    df.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/TCGA/ASE_summary/" + www + "_AF_vs_position_translocated_autosome.csv")

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

        if temp_chr == "chrX":
            xranges, colors = [(0,4400000)], [(.7, .7, .7)]
        elif temp_chr == "chr1":
            xranges, colors = [(0,2300000)], [(.7, .7, .7)]
        elif temp_chr == "chr17":
            xranges, colors = [(0,3400000)], [(.7, .7, .7)]
        elif temp_chr == "chr19":
            xranges, colors = [(0,6900000)], [(.7, .7, .7)]
        elif temp_chr == "chr22":
            xranges, colors = [(0,4300000)], [(.7, .7, .7)]

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

    fig, axs = plt.subplots(4, 3, sharey='row', gridspec_kw={'height_ratios': [6, 12, 3, 1], 'width_ratios': [1, 3, 1]})
    d = {}
    yticks = []
    yticklabels = []

    # ideogram.txt downloaded from UCSC's table browser
    for xranges, yrange, colors, label in ideograms('/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/cytoBandIdeo_' + temp_chr + '.txt'):
        print(xranges)
        coll = BrokenBarHCollection(xranges, yrange, facecolors=colors, lw=0)
        axs[3,1].add_collection(coll)
        center = yrange[0] + yrange[1]/2.
        yticks.append(center)
        yticklabels.append(label)
        d[label] = xranges

    axs[3,1].axis('tight')
    axs[3,1].set_yticks([])
    axs[3,1].set_yticklabels([])
    axs[3,1].set_xticks([])

    axs[3,1].set_xlim([-1000000, temp_length+1000000])

    #Plot copy number seg
    axs[2,1].set_xlim([-1000000, temp_length+1000000])
    axs[2,1].set_ylim([-0.01, 6.01]) 
    axs[2,1].set_yticks([0, 1, 2, 3, 4, 5, 6]) 
    axs[2,1].set_yticklabels([0, 1, 2, 3, 4, 5, 6]) 

    if "Normal" not in value:
        sub_cn_df = df_abs_cn[df_abs_cn['Sample'].str.contains(temp_id+"-01")]
    else:
        sub_cn_df = df_abs_cn[df_abs_cn['Sample'].str.contains(temp_id+"-1")]
    sub_cn_df = sub_cn_df[sub_cn_df['Chromosome'].isin((temp_chr.split("chr")[-1],int(temp_chr.split("chr")[-1])))]
    start_vals = sub_cn_df['Start'].tolist()
    end_vals = sub_cn_df['End'].tolist()
    cn_seg = [round((2**x)*2) for x in sub_cn_df['Segment_Mean'].tolist()]
    
    if "CJ-5681" in value and "Normal" not in value:
        sub_cn_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/CNIs/TCGA-CJ-5681_Tumor_cluster2.segs.txt", sep="\t")
        sub_cn_df = sub_cn_df[sub_cn_df['Chromosome'].isin((temp_chr.split("chr")[-1],int(temp_chr.split("chr")[-1])))]
        start_vals = sub_cn_df['Start_Position.bp.'].tolist()
        end_vals = sub_cn_df['End_Position.bp.'].tolist()
        cn_seg = sub_cn_df['Corrected_Copy_Number'].tolist()
    elif "BP-4756" in value and "Normal" not in value:
        sub_cn_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/CNIs/TCGA-BP-4756_Tumor_cluster1.segs.txt", sep="\t")
        sub_cn_df = sub_cn_df[sub_cn_df['Chromosome'].isin((temp_chr.split("chr")[-1],int(temp_chr.split("chr")[-1])))]
        start_vals = sub_cn_df['Start_Position.bp.'].tolist()
        end_vals = sub_cn_df['End_Position.bp.'].tolist()
        cn_seg = sub_cn_df['Corrected_Copy_Number'].tolist()

    """
    #Plot SNP level CN plots for these samples
    if ("G7-7501" in value or "B0-5705" in value or "T7-A92I" in value) and "Normal" not in value:
        sub_cn_df = pd.read_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/TCGA/representative_sample_CN/" + value + "_Tumor_cluster1.titan.txt", sep="\t")
        #sub_cn_df = sub_cn_df[sub_cn_df['Depth']>=50]
        sub_cn_df['logR_alleleA_product'] = sub_cn_df['Corrected_Ratio']*sub_cn_df['logR_Copy_Number']
        sub_cn_df['logR_alleleB_product'] = sub_cn_df['logR_Copy_Number']-sub_cn_df['logR_alleleA_product']
        start_vals = sub_cn_df['Position'].tolist()
        end_vals = sub_cn_df['Position'].tolist()
        alleleA_CN = sub_cn_df['logR_alleleA_product'].tolist()
        alleleB_CN = sub_cn_df['logR_alleleB_product'].tolist()
        qx=0
        while qx<len(start_vals):
            axs[2,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(alleleA_CN[qx]), float(alleleA_CN[qx])), marker=None, linestyle="-", lw=3, c="black")
            axs[2,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(alleleB_CN[qx]), float(alleleB_CN[qx])), marker=None, linestyle="-", lw=3, c="orange")
            qx=qx+1
    else:
    """

    qx=0
    while qx<len(start_vals):
        axs[2,1].plot((float(start_vals[qx]), float(end_vals[qx])), (float(cn_seg[qx]), float(cn_seg[qx])), marker=None, linestyle="-", lw=3, c="black")
        qx=qx+1

    #Plot minor allele fraction vs. position

    axs[1,1].scatter(xs,ys, c="black", alpha=1, s=2)
    axs[1,1].set_facecolor("white")
    axs[1,1].grid(False)
    axs[1,1].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4) 
    axs[1,1].spines['bottom'].set_color('0')
    axs[1,1].spines['left'].set_color('0')

    #Plot CDF of MAF (expressed SNVs + biallelic SNVs)
    count, bins_count = np.histogram(xs, bins=5000)
    pdf = count / sum(count)
    cdf = np.cumsum(pdf)
    axs[0,1].plot(bins_count[1:], cdf, c="black")
    axs[0,1].set_xlim([-1000000, temp_length+1000000])
    axs[0,1].set_ylim([-0.01, 1.01])

    x_max = bins_count[np.argmax(cdf >= 1)]
    x_min = bins_count[np.argmax(cdf <= 0)]
    cdf_min = np.min(cdf)
    print(x_max)
    #axs[0,1].plot((max(xs), temp_length), (1, 1), linewidth=1.4, c="black", marker=None, dash_capstyle="butt")
    #axs[0,1].plot((0, min(xs)-500000), (0, cdf_min), linewidth=1.4, c="black", marker=None, dash_capstyle="butt")

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
    axs[0,1].plot(bins_count[1:], cdf, c="red")
    axs[0,1].set_xlim([-1000000, temp_length+1000000])
    axs[0,1].set_ylim([-0.01, 1.01])

    x_max = bins_count[np.argmax(cdf >= 1)]
    x_min = bins_count[np.argmax(cdf <= 0)]
    cdf_min = np.min(cdf)
    print(cdf_min)
    #axs[0,1].plot((max(sub_xs), temp_length), (1, 1), linewidth=1.4, c="red", marker=None, dash_capstyle="butt")
    #axs[0,1].plot((0, min(sub_xs)-500000), (0, cdf_min), linewidth=1.4, c="red", marker=None, dash_capstyle="butt")

    patch1 = mpatches.Patch(color="red", label="Biallelic\nSNVs (RNA)")
    patch2 = mpatches.Patch(color="black", label="Expressed\nSNVs")

    legend = axs[0,1].legend(handles=[patch1, patch2], bbox_to_anchor=(-0.25, 0.84))

    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')

    axs[0,1].set_facecolor("white")
    axs[0,1].grid(False)
    axs[0,1].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4) 
    axs[0,1].tick_params(axis='y', which='major', left=True, labelleft=True, labelsize=10, size=4) 
    axs[0,1].spines['bottom'].set_color('0')
    axs[0,1].spines['left'].set_color('0')

    #Plot MAF distributions left and right of break
    distal = df[(df['position']<temp_loc)]
    proximal = df[(df['position']>=temp_loc)]
    distal_vals = [x for x in distal['MAF'].tolist() if x==x]
    prox_vals = [x for x in proximal['MAF'].tolist() if x==x]

    try:
        wz, p = scipy.stats.mannwhitneyu(prox_vals,distal_vals, alternative="two-sided")
    except:
        p = 1
    
    KS_pval.append(p)
    distal_vals_total.append(np.mean(distal_vals))
    prox_vals_total.append(np.mean(prox_vals))

    axs[3,1].set_title('P = %s' % ("{:.2E}".format(p)), fontdict={'fontsize': 18}, y=-1.5)

    axs[1,1].set_ylim([-0.004, 0.504])
    axs[1,0].set_ylim([-0.004, 0.504])
    axs[1,2].set_ylim([-0.004, 0.504])

    k1 = sns.kdeplot(np.array(distal_vals), weights=np.ones(len(distal_vals)) / len(distal_vals), bw=0.5, ax=axs[1,0], vertical=True, color="red", lw=1.5, fill="red")
    k2 = sns.kdeplot(np.array(prox_vals), weights=np.ones(len(prox_vals)) / len(prox_vals), bw=0.5, ax=axs[1,2], vertical=True, color="blue", lw=1.5, fill="blue")

    lim1 = axs[1,0].get_xlim()[1]
    lim2 = axs[1,0].get_xlim()[0]
    lim3 = axs[1,2].get_xlim()[1]
    lim4 = axs[1,2].get_xlim()[0]
    max_lim = max(lim1,lim2,lim3,lim4)
    axs[1,0].set_xlim([0,max_lim])
    axs[1,2].set_xlim([0,max_lim]) 
 
    axs[1,0].invert_xaxis()

    axs[1,1].axvspan(0, temp_loc, alpha=0.15, color='red')
    axs[1,1].axvspan(temp_loc, temp_length, alpha=0.15, color='blue')

    #Label axes
    axs[0,1].set_ylabel("Cumulative Probability")
    axs[0,1].set_xlabel("hg38 %s Position (bp)" % (temp_chr), labelpad=5)
    axs[1,1].set_ylabel("RNA Minor Allele Fraction", labelpad=5)
    axs[1,1].set_xlabel("hg38 %s Position (bp)" % (temp_chr), labelpad=5)
    axs[1,0].set_xlabel("Density\nLeft of %s" % (temp_partner))
    axs[1,2].set_xlabel("Density\nRight of %s" % (temp_partner))
    axs[2,1].set_ylabel("Copy\nNumber")
    axs[2,1].set_xlabel("hg38 %s Position (bp)" % (temp_chr), labelpad=5)
    axs[1,1].yaxis.tick_right()
    axs[1,1].set_xlim([-1000000, temp_length+1000000])

    #Formatting
    axs[1,0].set_facecolor("white")
    axs[1,0].grid(False)
    axs[1,0].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4)
    axs[1,0].tick_params(axis='y', which='major', right=True, labelsize=10, size=4, labelright=False, labelleft=False)
    axs[1,0].spines['bottom'].set_color('0')
    axs[1,0].spines['right'].set_color('0')

    axs[1,2].set_facecolor("white")
    axs[1,2].grid(False)
    axs[1,2].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4)
    axs[1,2].tick_params(axis='y', which='major', left=True, labelsize=10, size=4)
    axs[1,2].spines['bottom'].set_color('0')
    axs[1,2].spines['left'].set_color('0')

    axs[2,1].set_facecolor("white")
    axs[2,1].grid(False)
    axs[2,1].tick_params(axis='x', which='major', bottom=True, labelsize=10, size=4)
    axs[2,1].tick_params(axis='y', which='major', left=True, labelsize=10, size=4, labelleft=True)
    axs[2,1].spines['bottom'].set_color('0')
    axs[2,1].spines['left'].set_color('0')

    #axs[0,1].spines['bottom'].set_color('0')
    #axs[0,1].spines['left'].set_color('0')

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

    axs[3,2].set_facecolor("white")
    axs[3,2].grid(False)

    axs[3,0].set_facecolor("white")
    axs[3,0].grid(False)

    axs[3,2].get_xaxis().set_visible(False)
    axs[3,2].get_yaxis().set_visible(False)
    axs[3,0].get_xaxis().set_visible(False)
    axs[3,0].get_yaxis().set_visible(False)

    axs[0,1].axvline(temp_loc, c="black", lw=0.5, ls="--")
    axs[1,1].axvline(temp_loc, c="black", lw=0.5, ls="--")
    axs[2,1].axvline(temp_loc, c="black", lw=0.5, ls="--")

    fig.tight_layout()
    plt.subplots_adjust(wspace=0.55)
    plt.subplots_adjust(hspace=0.44)

    axs[1,1].tick_params(axis='y', which='major', left=True, right=True, labelsize=10, size=4, labelleft=True, labelright=True)

    fig.savefig("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/TCGA/unpaired_RNAseq_plots/fusionPARTNER_" + www + "_AF_vs_position.pdf", dpi=dpi_set, bbox_inches = 'tight')

df_out = pd.DataFrame([vals, distal_vals_total, prox_vals_total, KS_pval]).T
df_out.columns = ['Sample', 'Left_Mean_MAF', 'Right_Mean_MAF', 'pval_KS']
df_out.to_csv("/Users/ananthansadagopan/Documents/ViswanathanLab/tRCC_X/TCGA/ASE_summary/TCGA_translocated_autosome_KS_test.csv")
