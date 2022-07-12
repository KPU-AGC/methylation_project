#!/bin/bash
#Set correct directories
data=../data
results=../results
coverage=$data/coverage
bspcr_bed=$data/BSPCR_regions.bed

for file in $coverage/*.cov; do 
    echo $file
    python ./process-bismark-coverage.py $file $bspcr_bed -o $results
done