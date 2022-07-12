#!/bin/bash
#Set correct directories
data=../data
results=../results
trimmed=$results/trimmed
mapped=$results/mapped
fastqc=$results/fastqc
conda_env=~/miniconda3/envs/ngs/share/trimmomatic-0.39-2
filter_incomplete=$results/incomplete_conversion
index=$results/indexed
bisulfite_genome=/data/ref/genomes/arsucd1.2/slugger-masked
methylation_extraction=$results/methylation_extraction
#Check if directories exist, if not create
[ ! -d $fastqc ] && mkdir -p $fastqc
[ ! -d $trimmed ] && mkdir -p $trimmed
[ ! -d $mapped ] && mkdir -p $mapped
[ ! -d $filter_incomplete ] && mkdir -p $filter_incomplete
[ ! -d $index ] && mkdir -p $index

#STEP 1: TRIMMING AND QC
#Get fastqc results
fastqc -t 16 -outdir $results/fastqc $data/*
#Get prefixes, which correspond to each ngs run
ls $data | grep fastq | sed s/_R1_001.fastq//g | sed s/_R2_001.fastq//g | sort | uniq > $data/ngs_runlist.txt
#Trimming
while IFS= read -r prefix; do
       trimmomatic PE -threads 16 -phred33 -trimlog $trimmed/"$prefix"_trimlog.log $data/"$prefix"_R1_001.fastq $data/"$prefix"_R2_001.fastq $trimmed/"$prefix"_1.fastq $trimmed/"$prefix"_un1.fastq $trimmed/"$prefix"_2.fastq $trimmed/"$prefix"_un2.fastq ILLUMINACLIP:$conda_env/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:15
done < $data/ngs_runlist.txt
#Get fastqc results on trimmed sequences
fastqc -t 16 -outdir $results/fastqc $trimmed/*.fastq

#STEP 2: MAPPING
while IFS= read -r prefix; do
        bismark --genome $bisulfite_genome -1 $trimmed/"$prefix"_1.fastq -2 $trimmed/"$prefix"_2.fastq --non_directional -o $mapped
        filter_non_conversion --paired $mapped/"$prefix"_1_bismark_bt2_pe.bam
done < $data/ngs_runlist.txt
#Move files to correct directory
mv -t $filter_incomplete $mapped/*nonCG_filtered* $mapped/*.non-conversion_filtering.txt $mapped/*.nonCG_removed_seqs.bam

#STEP 3: INDEXING FOR VIEWING IN IGV
while IFS= read -r prefix; do
    samtools sort $filter_incomplete/"$prefix"_1_bismark_bt2_pe.nonCG_filtered.bam -o $index/"$prefix".sorted.bam
    samtools index $index/"$prefix".sorted.bam
done < $data/ngs_runlist.txt

#STEP 4: METHYLATION EXTRACTION
while IFS= read -r prefix; do
    bismark_methylation_extractor --gzip --bedGraph --cutoff 50 -p --parallel 4 -o $methylation_extraction $filter_incomplete/"$prefix"_1_bismark_bt2_pe.nonCG_filtered.bam 
done < $data/ngs_runlist.txt
