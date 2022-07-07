#!/bin/bash
#Set correct directories
data=../data
results=../results
filter_incomplete=$results/incomplete_conversion
methylation_extraction=$results/methylation_extraction
#Check if directories exist, if not create
[ ! -d $methylation_extraction ] && mkdir -p $methylation_extraction

#Methylation extraction
while read prefix
do
    bismark_methylation_extractor --gzip --bedGraph --cutoff 30 -p --parallel 4 -o $methylation_extraction $filter_incomplete/"$prefix"_1_bismark_bt2_pe.nonCG_filtered.bam 
done < $data/ngs_runlist.txt

mv -t $methylation_extraction $filter_incomplete/*.gz ß