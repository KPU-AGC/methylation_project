#!/bin/bash
#Set correct directories
#Main directories
data=../data
results=../results

#Trimming and QC
trimmed=$results/trimmed_bam
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

methylation_extraction=$bismark/methylation_extraction
all_extracted=$methylation_extraction/nonCG_filtered
filtered_extracted=$methylation_extraction/CG_filtered

#Default: /data/ref/genomes/arsucd1.2/slugger-masked
bisulfite_genome=$1

#Check if directories exist, if not create
[ ! -d $fastqc ] && mkdir -p $fastqc
[ ! -d $pre_trim_fastqc ] && mkdir -p $pre_trim_fastqc
[ ! -d $post_trim_fastqc ] && mkdir -p $post_trim_fastqc
[ ! -d $trimmed ] && mkdir -p $trimmed
[ ! -d $mapped ] && mkdir -p $mapped
[ ! -d $filter_incomplete ] && mkdir -p $filter_incomplete
[ ! -d $index ] && mkdir -p $index
[ ! -d $all_index ] && mkdir -p $all_index
[ ! -d $filtered_index ] && mkdir -p $filtered_index
[ ! -d $methylation_extraction ] && mkdir -p $methylation_extraction
[ ! -d $all_extracted ] && mkdir -p $all_extracted
[ ! -d $filtered_extracted ] && mkdir -p $filtered_extracted

#STEP 1: TRIMMING AND QC
#Get fastqc results for pre and post-trim
fastqc -t 16 -outdir $pre_trim_fastqc $data/*
#Pipeline for getting the sample name prefixes
#ls $data | grep fastq | sed s/_R1_001.fastq//g | sed s/_R2_001.fastq//g | sed s/_.*//g | sort | uniq > $data/ngs_samplelist.txt
ls $data | grep fastq | sed s/_R1_001.fastq//g | sed s/_R2_001.fastq//g | sort | uniq > $data/ngs_samplelist.txt
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
#Assumes that the bowtie2 mapper is used
while IFS= read -r prefix; do
        bismark --genome $bisulfite_genome -1 $trimmed/"$prefix"_1.fastq -2 $trimmed/"$prefix"_2.fastq --non_directional -o $mapped
        filter_non_conversion --paired $mapped/"$prefix"_1_bismark_bt2_pe.bam
done < $data/ngs_samplelist.txt
mv -t $filter_incomplete $mapped/*nonCG_filtered* $mapped/*.non-conversion_filtering.txt $mapped/*.nonCG_removed_seqs.bam

#STEP 3: INDEXING FOR VIEWING IN IGV
#Both the filtered and non-filtered BAMs will be used
while IFS= read -r prefix; do
    samtools sort $filter_incomplete/"$prefix"_1_bismark_bt2_pe.nonCG_filtered.bam -o $all_index/"$prefix".nonCG_filtered.sorted.bam
    samtools index $all_index/"$prefix".nonCG_filtered.sorted.bam
done < $data/ngs_samplelist.txt

while IFS= read -r prefix; do
    samtools sort $mapped/"$prefix"_1_bismark_bt2_pe.bam -o $filtered_index/"$prefix".sorted.bam
    samtools index $filtered_index/"$prefix".sorted.bam
done < $data/ngs_samplelist.txt

#STEP 4: METHYLATION EXTRACTION
#Both filtered and non-filtered BAMs will be used
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

#STEP 5: GENERATE REPORTS
