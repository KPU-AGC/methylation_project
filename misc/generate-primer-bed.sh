#!/bin/bash
#generate-primer-bed.sh : script to generate primer bed file for analysis scripts
#Additionally, it will generate fasta files for each of the primer regions of the bed file
#for troubleshooting purposes.

#Set correct directories
#CHANGE genome and output_name TO FIT WHAT YOU NEED
data=../data
results=../results
fasta=$results/fasta
#input_name=$data/bspcr-primers.tab
input_name=$1
#genome=$data/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa
genome=$2
bed_output=$results/BSPCR-primers.bed
fasta_output=$results/BSPCR-all-primers.fasta

[ ! -d $fasta ] && mkdir -p $fasta

bedtools getfasta -fi $genome -bed $input_name -bedOut > $bed_output
bedtools getfasta -fi $genome -bed $input_name -name > $fasta_output

python ./split-fasta.py $fasta_output --output $fasta