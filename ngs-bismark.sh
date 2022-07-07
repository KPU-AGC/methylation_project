#!/bin/bash
#Set correct directories
data=../data
results=../results
trimmed=$results/trimmed
mapped=$results/mapped
filter_incomplete=$results/incomplete_conversion
bisulfite_genome=/data/ref/genomes/arsucd1.2/slugger-masked
#Check if directories exist, if not create
[ ! -d $mapped ] && mkdir -p $mapped
[ ! -d $filter_incomplete ] && mkdir -p $filter_incomplete

#Read Mapping
while read prefix
do
        bismark --genome $bisulfite_genome -1 $trimmed/"$prefix"_1.fastq -2 $trimmed/"$prefix"_2.fastq --non_directional -o $mapped
        filter_non_conversion --paired $mapped/"$prefix"_1_bismark_bt2_pe.bam
        #samtools sort $mapped/"$prefix"_1_bismark_bt2_pe.bam -o $mapped/"$prefix"_1_bismark_bt2_pe-sorted.bam
        #samtools index $mapped/"$prefix"_1_bismark_bt2_pe-sorted.bam
done < $data/ngs_runlist.txt
#Move files to correct directory
mv -t $filter_incomplete $mapped/*nonCG_filtered* $mapped/*.non-conversion_filtering.txt $mapped/*.nonCG_removed_seqs.bam
#mv *nonCG_filtered* $filter_incomplete
#mv *.non-conversion_filtering.txt $filter_incomplete