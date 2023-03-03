#!/bin/bash
#get-fasta.sh : scripts to extract individual fasta sequences from a genome at specified genomic coordinates. 

#Set correct directories
#CHANGE genome and output_name TO FIT WHAT YOU NEED
data=../data
results=../results
fasta=$results/fasta
input_name=$data/bspcr-primers.tab
genome=$data/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa
bed_output=$results/BSPCR-primers.bed
fasta_output=$results/BSPCR-all-primers.fasta

[ ! -d $fasta ] && mkdir -p $fasta

bedtools getfasta -fi $genome -bed $input_name -bedOut > $bed_output
bedtools getfasta -fi $genome -bed $input_name -name > $fasta_output

python ./split-fasta.py $fasta_output --output $fasta