#!/bin/bash
#ngs-process.sh : script to run pipeline to process FASTQ files from targeted bisulfite amplicon sequencing 
#into coverage and bedgraph files for downstream analysis.
#Need to pass path to the reference genome for the analysis to be performed.
#----------------------------------------------------------------------------------------------------------
#Set correct directories - directories will be made at each step they are used.
#NOTE: directories relative to script being location in bin folder. 
#See AGC document on data analysis directory tree (link)
#------------------------------------------------------------------------------

#STEP 0: SETTING UP DIRECTORIES
#------------------------------
data=../data
results=../results
#Trimming and QC
trimmed=$results/trimmed_fastq
fastqc=$results/fastqc
pre_trim_fastqc=$fastqc/pre_trim
post_trim_fastqc=$fastqc/post_trim
#conda_env needs to be changed to the correct location for the environment you're using
conda_env=~/miniconda3/envs/ngs/share/trimmomatic-0.39-2
#Bismark
bismark=$results/bismark
mapped=$bismark/mapped_bam
filter_incomplete=$bismark/incomplete_conv_bam
index=$bismark/indexed_bam
all_index=$index/nonCG_filtered
filtered_index=$index/CG_filtered
#Bismark methylation extraction
methylation_extraction=$bismark/methylation_extraction
all_extracted=$methylation_extraction/nonCG_filtered
filtered_extracted=$methylation_extraction/CG_filtered

#Default: /data/ref/genomes/arsucd1.2/slugger-masked
bisulfite_genome=$1

#STEP 1: TRIMMING AND QC
#-----------------------
#Generating the directories
[ ! -d $fastqc ] && mkdir -p $fastqc
[ ! -d $pre_trim_fastqc ] && mkdir -p $pre_trim_fastqc
[ ! -d $post_trim_fastqc ] && mkdir -p $post_trim_fastqc

#Get fastqc results for pre-trim files
fastqc -t 16 -outdir $pre_trim_fastqc $data/*

#Pipeline for getting the sample name prefixes
#ls $data | grep fastq | sed s/_R1_001.fastq//g | sed s/_R2_001.fastq//g | sed s/_.*//g | sort | uniq > $data/ngs_samplelist.txt
ls $data | grep fastq | sed s/_R1_001.fastq//g | sed s/_R2_001.fastq//g | sort | uniq > $data/ngs_samplelist.txt

#Generating the directories for trimmed output
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

#STEP 2: MAPPING
#---------------
#Assumes that the bowtie2 mapper is used

#Generating directory for mapped output
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

[ ! -d $index ] && mkdir -p $index
[ ! -d $all_index ] && mkdir -p $all_index
[ ! -d $filtered_index ] && mkdir -p $filtered_index
#Non-filtered - BAM sequences need to be sorted before being indexed
while IFS= read -r prefix; do
    samtools sort $filter_incomplete/"$prefix"_1_bismark_bt2_pe.nonCG_filtered.bam -o $all_index/"$prefix".nonCG_filtered.sorted.bam
    samtools index $all_index/"$prefix".nonCG_filtered.sorted.bam
done < $data/ngs_samplelist.txt
#Filtered
while IFS= read -r prefix; do
    samtools sort $mapped/"$prefix"_1_bismark_bt2_pe.bam -o $filtered_index/"$prefix".sorted.bam
    samtools index $filtered_index/"$prefix".sorted.bam
done < $data/ngs_samplelist.txt

#STEP 4: METHYLATION EXTRACTION
#------------------------------
#Perform for the filtered and unfiltered BAM files
[ ! -d $methylation_extraction ] && mkdir -p $methylation_extraction
[ ! -d $all_extracted ] && mkdir -p $all_extracted
[ ! -d $filtered_extracted ] && mkdir -p $filtered_extracted

#Non-filtered
while IFS= read -r prefix; do
    bismark_methylation_extractor \
        --gzip \
        --bedGraph \
        --cutoff 50 \
        -p \
        --parallel 4 \
        -o $all_extracted \
        $filter_incomplete/"$prefix"_1_bismark_bt2_pe.nonCG_filtered.bam 
done < $data/ngs_samplelist.txt
#Filtered
while IFS= read -r prefix; do
    bismark_methylation_extractor \
        --gzip \
        --bedGraph \
        --cutoff 50 \
        -p \
        --parallel 4 \
        -o $filtered_extracted \
        $mapped/"$prefix"_1_bismark_bt2_pe.bam
done < $data/ngs_samplelist.txt