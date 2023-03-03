#!/bin/bash
#ngs-qc.sh : script to run fastqc on everything.

#Set correct directories
data=../data
results=../results
trimmed=$results/trimmed
fastqc=$results/fastqc
conda_env=~/miniconda3/envs/ngs-qc/share/trimmomatic-0.39-2
#Check if directories exist, if not create
[ ! -d $fastqc ] && mkdir -p $fastqc
[ ! -d $trimmed ] && mkdir -p $trimmed

#Get fastqc results
fastqc -t 16 -outdir $results/fastqc $data/*
#Get prefixes, which correspond to each ngs run
ls $data | grep fastq | sed s/_R1_001.fastq//g | sed s/_R2_001.fastq//g | sort | uniq > $data/ngs_runlist.txt
#Trimming
while read prefix
do
       trimmomatic PE -threads 16 -phred33 -trimlog $trimmed/"$prefix"_trimlog.log $data/"$prefix"_R1_001.fastq $data/"$prefix"_R2_001.fastq $trimmed/"$prefix"_1.fastq $trimmed/"$prefix"_un1.fastq $trimmed/"$prefix"_2.fastq $trimmed/"$prefix"_un2.fastq ILLUMINACLIP:$conda_env/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:4:15
done < $data/ngs_runlist.txt
#Get fastqc results on trimmed sequences
fastqc -t 16 -outdir $results/fastqc $trimmed/*.fastq