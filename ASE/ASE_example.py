#qsub -t 1:1 ASE_submit.qsub
import glob, os
from os import listdir
from os.path import isfile, join
import pandas as pd
import gzip

#Samples to iterate over, have DNA/RNA downloaded for these samples
vals = ['B_TRCC_18_Tumor']

def generate_vals(path1):

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

    names, header_temp = get_vcf_names(path1)
    df = pd.read_csv(path1, comment='#', delim_whitespace=True, header=None, names=names, compression="gzip")

    df['l1'] = df.REF.str.len()
    df['l2'] = df.ALT.str.len()

    df = df[df['l1']==1]
    df = df[df['l2']==1]

    del df['l1']
    del df['l2']

    file2 = path1[:-7] + "_filtered_temp.vcf"
    file3 = path1[:-7] + "_temp.txt"
    file4 = path1[:-7] + "_filtered.vcf"
    df.to_csv(file2, index=False, sep="\t")

    with open(file3, "w+") as file1:
        for q in header_temp:
            file1.write(q)

    os.system(("cat %s %s > %s") % (file3,file2,file4))

task_id = int(os.getenv("SGE_TASK_ID"))
param1 = vals[task_id-1]

def generate(param1):

    #RNA processing
    os.system("samtools sort -n /mnt/scratch/DFCI_tRCCs/%s/%s_RNA.bam -o /mnt/scratch/DFCI_tRCCs/%s/%s_RNA_sorted.bam" % (param1,param1,param1,param1))
    os.system("bedtools bamtofastq -i /mnt/scratch/DFCI_tRCCs/%s/%s_RNA_sorted.bam -fq /mnt/scratch/DFCI_tRCCs/%s/%s_RNA_1.fq -fq2 /mnt/scratch/DFCI_tRCCs/%s/%s_RNA_2.fq" % (param1,param1,param1,param1,param1,param1))
    os.system("tophat --mate-inner-dist 100 --output-dir /mnt/scratch/DFCI_tRCCs/%s/ /mnt/storage/labs/sviswanathan/GATK_ref/bowtie_index /mnt/scratch/DFCI_tRCCs/%s/%s_RNA_1.fq /mnt/scratch/DFCI_tRCCs/%s/%s_RNA_2.fq" % (param1,param1,param1,param1,param1))
    os.system("samtools sort /mnt/scratch/DFCI_tRCCs/%s/accepted_hits.bam -o /mnt/scratch/DFCI_tRCCs/%s/accepted_hits_sorted.bam" % (param1,param1))
    os.system("java -jar /home/au409/picard/build/libs/picard.jar AddOrReplaceReadGroups I=/mnt/scratch/DFCI_tRCCs/%s/accepted_hits_sorted.bam O=/mnt/scratch/DFCI_tRCCs/%s/accepted_hits_sorted_picard.bam RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=CELL VALIDATION_STRINGENCY=LENIENT" % (param1,param1))
    os.system("samtools index /mnt/scratch/DFCI_tRCCs/%s/accepted_hits_sorted_picard.bam" % (param1))

    #DNA processing
    os.system("samtools sort -n /mnt/scratch/DFCI_tRCCs/%s/%s_DNA.bam -o /mnt/scratch/DFCI_tRCCs/%s/%s_DNA_sorted.bam" % (param1,param1,param1,param1))
    os.system("bedtools bamtofastq -i /mnt/scratch/DFCI_tRCCs/%s/%s_DNA_sorted.bam -fq /mnt/scratch/DFCI_tRCCs/%s/%s_DNA_1.fq -fq2 /mnt/scratch/DFCI_tRCCs/%s/%s_DNA_2.fq" % (param1,param1,param1,param1,param1,param1))
    os.system("/mnt/storage/apps/BWA/bwa-0.7.17/bwa mem /mnt/storage/labs/sviswanathan/GATK_ref/Homo_sapiens_assembly38.fasta /mnt/scratch/DFCI_tRCCs/%s/%s_DNA_1.fq /mnt/scratch/DFCI_tRCCs/%s/%s_DNA_2.fq -o /mnt/scratch/DFCI_tRCCs/%s/%s_hg38.sam" % (param1,param1,param1,param1,param1,param1))
    os.system("samtools view -S -b /mnt/scratch/DFCI_tRCCs/%s/%s_hg38.sam > /mnt/scratch/DFCI_tRCCs/%s/%s_hg38.bam" % (param1,param1,param1,param1))
    os.system("samtools sort /mnt/scratch/DFCI_tRCCs/%s/%s_hg38.bam -o /mnt/scratch/DFCI_tRCCs/%s/%s_hg38_sorted.bam" % (param1,param1,param1,param1))
    os.system("java -jar /home/au409/picard/build/libs/picard.jar AddOrReplaceReadGroups I=/mnt/scratch/DFCI_tRCCs/%s/%s_hg38_sorted.bam O=/mnt/scratch/DFCI_tRCCs/%s/%s_hg38_sorted_picard.bam RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=CELL VALIDATION_STRINGENCY=LENIENT" % (param1,param1,param1,param1))
    os.system("samtools index /mnt/scratch/DFCI_tRCCs/%s/%s_hg38_sorted_picard.bam" % (param1,param1))
    os.system("java -Djava.io.tmpdir=/mnt/scratch -Xmx8000m -jar /mnt/storage/apps/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar HaplotypeCaller -R /mnt/storage/labs/sviswanathan/GATK_ref/Homo_sapiens_assembly38.fasta -I /mnt/scratch/DFCI_tRCCs/%s/%s_hg38_sorted_picard.bam -O /mnt/scratch/DFCI_tRCCs/%s/%s_hg38_sorted_picard.g.vcf.gz -ERC GVCF" % (param1,param1,param1,param1))
    os.system("java -Djava.io.tmpdir=/mnt/scratch -jar /mnt/storage/apps/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar GenotypeGVCFs -R /mnt/storage/labs/sviswanathan/GATK_ref/Homo_sapiens_assembly38.fasta -V /mnt/scratch/DFCI_tRCCs/%s/%s_hg38_sorted_picard.g.vcf.gz -O /mnt/scratch/DFCI_tRCCs/%s/%s_hg38_sorted_picard.vcf.gz" % (param1,param1,param1,param1))
    os.system("tabix -p vcf /mnt/scratch/DFCI_tRCCs/%s/%s_hg38_sorted_picard.vcf.gz" % (param1,param1))
    generate_vals("/mnt/scratch/DFCI_tRCCs/" + param1 + "/" + param1 + "_hg38_sorted_picard.vcf.gz")
    os.system("bgzip -c /mnt/scratch/DFCI_tRCCs/%s/%s_hg38_sorted_picard_filtered.vcf > /mnt/scratch/DFCI_tRCCs/%s/%s_hg38_sorted_picard_filtered.vcf.gz" % (param1,param1,param1,param1))
    os.system("tabix -p vcf /mnt/scratch/DFCI_tRCCs/%s/%s_hg38_sorted_picard_filtered.vcf.gz" % (param1,param1))
    
    #ASE
    os.system("java -Djava.io.tmpdir=/mnt/scratch/ -Xmx8000m -jar /mnt/storage/labs/sviswanathan/GATK_ref/cz_ref/libs/gatk.jar ASEReadCounter -R /mnt/storage/labs/sviswanathan/GATK_ref/Homo_sapiens_assembly38.fasta --read-filter PassesVendorQualityCheckReadFilter --read-filter HasReadGroupReadFilter --read-filter NotDuplicateReadFilter --read-filter MappingQualityAvailableReadFilter --read-filter NotSecondaryAlignmentReadFilter --read-filter MappingQualityReadFilter --minimum-mapping-quality 30 --read-filter OverclippedReadFilter --filter-too-short 25 --read-filter GoodCigarReadFilter --read-filter AmbiguousBaseReadFilter -V /mnt/scratch/DFCI_tRCCs/%s/%s_hg38_sorted_picard_filtered.vcf.gz --lenient --seconds-between-progress-updates 100 -L chrX -I /mnt/scratch/DFCI_tRCCs/%s/accepted_hits_sorted_picard.bam -O /mnt/scratch/DFCI_tRCCs/%s/RNA_paired_on_germline_hets_filtered_unphased_chrX_%s.log" % (param1,param1,param1,param1,param1))

generate(param1)

