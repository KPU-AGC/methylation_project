#!/bin/bash
#ngs-process.sh : script to run pipeline to process FASTQ files from targeted bisulfite amplicon sequencing 
#into coverage and bedgraph files for downstream analysis.
#PARAMETERS
#----------
#bisulfite_genome : path to folder containing the bisulfite-converted genome
#snp_file : path to SNP file necessary for SNPSplit to run. Lists one set of alleles that sequences will be separated into.
#----------------------------------------------------------------------------------------------------------
#Set correct directories - directories will be made at each step they are used.
#NOTE: directories relative to script being location in bin folder. 
#See AGC document on data analysis directory tree (link)
#------------------------------------------------------------------------------

#STEP 0: SETTING UP DIRECTORIES
#------------------------------
#IMPORTANT: conda_env needs to be changed to the correct location for the environment you're using
conda_env=~/miniconda3/envs/ngsv2/share/trimmomatic-0.39-2
data=../data
results=../results

#Bismark methylation extraction
methylation_extraction=$bismark/methylation_extraction
all_extracted=$methylation_extraction/all_sequences
filtered_extracted=$methylation_extraction/nonCG_filtered

#Default: /data/ref/genomes/arsucd1.2/slugger-masked
bisulfite_genome=$1
snp_file=$2

#STEP 1: TRIMMING AND QC
#-----------------------
#Generating the directories
fastqc=$results/fastqc
pre_trim_fastqc=$fastqc/pre_trim
post_trim_fastqc=$fastqc/post_trim
[ ! -d $fastqc ] && mkdir -p $fastqc
[ ! -d $pre_trim_fastqc ] && mkdir -p $pre_trim_fastqc
[ ! -d $post_trim_fastqc ] && mkdir -p $post_trim_fastqc

#Get fastqc results for pre-trim files
fastqc -t 16 -outdir $pre_trim_fastqc $data/*.fastq

#Pipeline for getting the sample name prefixes
#ls $data | grep fastq | sed s/_R1_001.fastq//g | sed s/_R2_001.fastq//g | sed s/_.*//g | sort | uniq > $data/ngs_samplelist.txt
ls $data | grep fastq | sed s/_R1_001.fastq//g | sed s/_R2_001.fastq//g | sort | uniq > $data/ngs_samplelist.txt

#Generating the directories for trimmed output
trimmed=$results/trimmed_fastq
[ ! -d $trimmed ] && mkdir -p $trimmed
#Trimming
while IFS= read -r prefix; do
       trimmomatic PE \
        -threads 16 \
        -phred33 \
        -trimlog $trimmed/"$prefix"_trimlog.log \
        $data/"$prefix"_R1_001.fastq \
        $data/"$prefix"_R2_001.fastq \
        $trimmed/"$prefix"_1.fastq \
        $trimmed/"$prefix"_un1.fastq \
        $trimmed/"$prefix"_2.fastq \
        $trimmed/"$prefix"_un2.fastq \
        ILLUMINACLIP:$conda_env/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:15
done < $data/ngs_samplelist.txt
fastqc -t 16 -outdir $post_trim_fastqc $trimmed/*.fastq

#Running multiqc
multiqc $pre_trim_fastqc --outdir $pre_trim_fastqc
multiqc $post_trim_fastqc --outdir $post_trim_fastqc

#STEP 2: MAPPING
#---------------
#Assumes that the bowtie2 mapper is used

#Generating directory for mapped output
bismark=$results/bismark
mapped=$bismark/mapped_bam
filter_incomplete=$bismark/incomplete_conv_bam
[ ! -d $mapped ] && mkdir -p $mapped
[ ! -d $filter_incomplete ] && mkdir -p $filter_incomplete

#Begin mapping - also perform non-converted filtering after mapping
while IFS= read -r prefix; do
        bismark --genome $bisulfite_genome -1 $trimmed/"$prefix"_1.fastq -2 $trimmed/"$prefix"_2.fastq --non_directional -o $mapped
        filter_non_conversion --paired $mapped/"$prefix"_1_bismark_bt2_pe.bam
done < $data/ngs_samplelist.txt
#Move results to the appropriate folder
mv -t $filter_incomplete $mapped/*nonCG_filtered* $mapped/*.non-conversion_filtering.txt $mapped/*.nonCG_removed_seqs.bam

#STEP 3: INDEXING FOR VIEWING IN IGV
#-----------------------------------
#Perform for the filtered and unfiltered BAM files

#Generating directories
index=$bismark/indexed_bam
all_index=$index/all_sequences
filtered_index=$index/nonCG_filtered
[ ! -d $index ] && mkdir -p $index
[ ! -d $all_index ] && mkdir -p $all_index
[ ! -d $filtered_index ] && mkdir -p $filtered_index

#Filtered - BAM sequences need to be sorted before being indexed
while IFS= read -r prefix; do
    samtools sort $filter_incomplete/"$prefix"_1_bismark_bt2_pe.nonCG_filtered.bam -o $filtered_index/"$prefix".nonCG_filtered.sorted.bam
    samtools index $filtered_index/"$prefix".nonCG_filtered.sorted.bam
done < $data/ngs_samplelist.txt

#Not filtered
while IFS= read -r prefix; do
    samtools sort $mapped/"$prefix"_1_bismark_bt2_pe.bam -o $all_index/"$prefix".sorted.bam
    samtools index $all_index/"$prefix".sorted.bam
done < $data/ngs_samplelist.txt

#STEP 4: SPLITTING BAM FILES BY PATERNAL AND MATERNAL ALLELES
#------------------------------------------------------------
#SNPsplit is used to separate BAM files using paternal genotypes.
#Maternal genotype is assumed to be the same as the reference genome.

#Generating directories
split=$bismark/split_bam
filtered_split=$split/nonCG_filtered
all_split=$split/all_sequences
[ ! -d $split ] && mkdir -p $split
[ ! -d $filtered_split ] && mkdir -p $filtered_split
[ ! -d $all_split ] && mkdir -p $all_split

#Filtered
while IFS= read -r prefix; do
    SNPsplit $filtered_index/"$prefix".nonCG_filtered.sorted.bam --paired --bisulfite --snp_file $snp_file -o $filtered_split
done < $data/ngs_samplelist.txt

#All sequences
while IFS= read -r prefix; do
    SNPsplit $all_index/"$prefix".sorted.bam --paired --bisulfite --snp_file $snp_file -o $all_split
done < $data/ngs_samplelist.txt

#STEP 5: METHYLATION EXTRACTION - NO SPLIT
#-----------------------------------------

#Generating directories
methylation_extraction=$bismark/methylation_extraction
all_extracted=$methylation_extraction/all_sequences
filtered_extracted=$methylation_extraction/nonCG_filtered

[ ! -d $methylation_extraction ] && mkdir -p $methylation_extraction
[ ! -d $all_extracted ] && mkdir -p $all_extracted
[ ! -d $filtered_extracted ] && mkdir -p $filtered_extracted

#Filtered
while IFS= read -r prefix; do
    bismark_methylation_extractor \
        --gzip \
        --bedGraph \
        --zero_based \
        -p \
        --parallel 4 \
        -o $filtered_extracted \
        $filter_incomplete/"$prefix"_1_bismark_bt2_pe.nonCG_filtered.bam 
done < $data/ngs_samplelist.txt
#Non-filtered
while IFS= read -r prefix; do
    bismark_methylation_extractor \
        --gzip \
        --bedGraph \
        --zero_based \
        -p \
        --parallel 4 \
        -o $all_extracted \
        $mapped/"$prefix"_1_bismark_bt2_pe.bam
done < $data/ngs_samplelist.txt

#STEP 6: METHYLATION EXTRACTION - 
snpsplit_extract=$bismark/snpsplit_methylation_extraction
all_snpsplit_extract=$snpsplit_extract/all_sequences
filtered_snpsplit_extract=$snpsplit_extract/nonCG_filtered
[ ! -d $snpsplit_extract ] && mkdir -p $snpsplit_extract
[ ! -d $all_snpsplit_extract ] && mkdir -p $all_snpsplit_extract
[ ! -d $filtered_snpsplit_extract ] && mkdir -p $filtered_snpsplit_extract

#Filtered
while IFS= read -r prefix; do
    #Genome1
    bismark_methylation_extractor \
        --gzip \
        --bedGraph \
        --zero_based \
        -p \
        --parallel 4 \
        -o $filtered_snpsplit_extract \
        $filtered_split/"$prefix".nonCG_filtered.sorted.genome1.bam
    #Genome2
    bismark_methylation_extractor \
        --gzip \
        --bedGraph \
        --zero_based \
        -p \
        --parallel 4 \
        -o $filtered_snpsplit_extract \
        $filtered_split/"$prefix".nonCG_filtered.sorted.genome2.bam
    #Unassigned
    bismark_methylation_extractor \
        --gzip \
        --bedGraph \
        --zero_based \
        -p \
        --parallel 4 \
        -o $filtered_snpsplit_extract \
        $filtered_split/"$prefix".nonCG_filtered.sorted.unassigned.bam
done < $data/ngs_samplelist.txt

#All sequences
while IFS= read -r prefix; do
    #Genome1
    bismark_methylation_extractor \
        --gzip \
        --bedGraph \
        --zero_based \
        -p \
        --parallel 4 \
        -o $all_snpsplit_extract \
        $all_split/"$prefix".sorted.genome1.bam
    #Genome2
    bismark_methylation_extractor \
        --gzip \
        --bedGraph \
        --zero_based \
        -p \
        --parallel 4 \
        -o $all_snpsplit_extract \
        $all_split/"$prefix".sorted.genome2.bam
    #Unassigned
    bismark_methylation_extractor \
        --gzip \
        --bedGraph \
        --zero_based \
        -p \
        --parallel 4 \
        -o $all_snpsplit_extract \
        $all_split/"$prefix".sorted.unassigned.bam
done < $data/ngs_samplelist.txt