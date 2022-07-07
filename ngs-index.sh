#!/bin/bash
#Set correct directories
data=../data
results=../results
filter_incomplete=$results/incomplete_conversion
index=$results/indexed
#Check if directories exist, if not create
[ ! -d $index ] && mkdir -p $index
while read prefix
do
    samtools sort $filter_incomplete/"$prefix"_1_bismark_bt2_pe.nonCG_filtered.bam -o $index/"$prefix".sorted.bam
    samtools index $index/"$prefix".sorted.bam
done < $data/ngs_runlist.txt