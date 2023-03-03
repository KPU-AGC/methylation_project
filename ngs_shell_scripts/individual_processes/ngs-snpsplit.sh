#!/bin/bash
#Set correct directories
data=../data
results=../results
filter_incomplete=$results/incomplete_conversion
snpsplit=$results/snpsplit
snp_file=$data/Slugger_variants.txt
#Check if directories exist, if not create
[ ! -d $snpsplit ] && mkdir -p $snpsplit

#Run snpsplit
while read prefix
do
    SNPsplit --paired --bisulfite --snp_file $snp_file $filter_incomplete/"$prefix"_1_bismark_bt2_pe.nonCG_filtered.bam
done < $data/ngs_runlist.txt

#Move to the correct folder
echo "Moving files..."
mv -t $snpsplit *.allele_flagged.bam *.genome1.bam *.genome2.bam *.SNPsplit_report.txt *.SNPsplit_report.yaml *.SNPsplit_sort.txt *.unassigned.bam
echo "Done!"