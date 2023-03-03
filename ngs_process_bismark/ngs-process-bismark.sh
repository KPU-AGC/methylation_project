#!/bin/bash
#Set correct directories
data=../data
results=../results
bspcr_bed=$data/BSPCR_regions.bed
sample_summary=$results/sample_summary
primer_summary=$results/primer_summary
detailed=$results/detailed

[ ! -d $sample_summary ] && mkdir -p $sample_summary
[ ! -d $primer_summary ] && mkdir -p $primer_summary
[ ! -d $detailed ] && mkdir -p $detailed

for file in $data/*.cov; do 
    echo $file
    python ./process-bismark-coverage.py $file $bspcr_bed -o $detailed
done

mv $detailed/*_summary.csv $sample_summary

python ./assess-by-primer.py $sample_summary --output $primer_summary