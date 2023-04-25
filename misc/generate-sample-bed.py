"""
get-sample-bed.py - Script to quickly generate a bed file as well as a SNPSplit input file. 

Date : 2023-04-24

Functions
---------
filter_variants : from the list of variants, filter based on selected criteria
generate_variant_bed : generate the BED file for the genome masking process
generate_snpsplit_file : generate the .tsv file used for SNPsplit

#Variant keys
# [0] : chromosome number
# [1] : nucleotide position on chromosome
# [2] : additional sequence ID field
# [3] : reference allele
# [4] : alt allele
# [5] : quality score
# [6] : filter value (LowQual, PASS)
# [7] : type (SNV, Deletion, Insertion)
# [8] : genotype (het., hom. ALT)
# [9] : base position in sequence file
# [10] : signal position

#chromosome number, nucleotide position-1, nucleotide position
#bed_format = (variant[0], variant[1]-1, variant[1])
#snp index, chromosome number, nucleotide position, 1, REF/ALT
#snp_format = (i, variant[0], variant[1], 1, f'{variant[3]}/{variant[4]}')             

"""

__author__ = 'Michael Ke'
__version__ = '1.1.0'

import argparse
import csv
import pathlib
import json

def parse_args(): 
    parser = argparse.ArgumentParser('Script to produce a BED file with variant positions')
    parser.add_argument(
        'json_path',
        action='store',
        type=pathlib.Path,
        help='Path to the json output of combine-variant-calls.py for a single sample.'
    )
    parser.add_argument(
        'target_sample',
        action='store',
        type=str,
        help='Key of target sample. This will be the prefix of output files.'
    )
    parser.add_argument(
        '--output',
        dest='output_path',
        action='store',
        default=None,
        type=pathlib.Path,
        help='Output path'
    )
    parser.add_argument(
        '--snp_only', 
        dest='snp_flag',
        action='store_true',
        help='Flag to only include SNPs'
    )
    parser.add_argument(
        '--no_het',
        dest='het_flag',
        action='store_true',
        help='Flag to only include homozygous variants'
    )
    args = parser.parse_args()

    if args.output_path is None: 
        args.output_path = args.json_path.parent
    
    return args

def filter_variants(data: list, snp_flag: bool, het_flag: bool) -> list: 
    """
    Filter variants according to flags.

    Parameters
    ----------
    data : list
        list of variant data
    snp_flag : bool
        If true, filter out non-SNP variants.
    het_flag : bool
        If true, only include homozygous variants.
    
    """

    filtered_variants = []

    for variant in data: 
        #Screen variants
        variant_include = True
        if snp_flag and not(variant[7] == 'SNV'): 
            variant_include = False
        elif het_flag and not(variant[8] == 'hom. ALT'): 
            variant_include = False
        else: 
            pass
        
        if variant_include: 
            filtered_variants.append(variant)

    return filtered_variants

def generate_variant_bed(data: list, name: str, output_path: pathlib.Path) -> None: 
    """Generate the BED file of variants for an individual."""

    bed_output_path = output_path.joinpath(f'{name}_variants.bed')

    #Extract relevant fields for BED file, and remove duplicates from 
    # FW and REV reads
    variant_set = set()
    for variant in data:
        variant_set.add(
            (
                variant[0], 
                variant[1]-1, 
                variant[1],
                f'{variant[3]}>{variant[4]}'
            )
        )
    
    with open(bed_output_path, 'w', newline='') as bed_file: 
        csv_writer = csv.writer(bed_file, delimiter='\t')
        csv_writer.writerows(variant_set)

def generate_snpsplit_file(data: list, name: str, output_path: pathlib.Path) -> None: 
    """Generate the file necessary for SNPsplit to work"""
    
    snp_output_path = output_path.joinpath(f'{name}_snpsplit.tsv')

    #Extract relevant fields for BED file, and remove duplicates from 
    # FW and REV reads
    variant_set = set()
    
    #The variants are arbitrarily indexed. They just need an ID name.
    for variant in data:
        variant_set.add(
            (
                variant[2].split('_')[1].split('-')[0],
                #Troubleshooting - show the exact primer set used
                #variant[2].split('_')[1],
                variant[0], 
                variant[1], 
                1, 
                f'{variant[3]}/{variant[4]}'
            )
        )
    
    with open(snp_output_path, 'w', newline='') as bed_file: 
        csv_writer = csv.writer(bed_file, delimiter='\t')
        csv_writer.writerows(variant_set) 

def main(): 
    args = parse_args()

    #Loading data
    json_file = open(args.json_path, 'r')
    json_data = json.load(json_file)
    json_file.close()

    #Filter variants, if necessary
    filtered_variants = filter_variants(json_data['variants'], args.snp_flag, args.het_flag)

    generate_variant_bed(filtered_variants, args.target_sample, args.output_path)

    generate_snpsplit_file(filtered_variants, args.target_sample, args.output_path)

    #Variant keys
    # [0] : chromosome number
    # [1] : nucleotide position on chromosome
    # [2] : additional sequence ID field
    # [3] : reference allele
    # [4] : alt allele
    # [5] : quality score
    # [6] : filter value (LowQual, PASS)
    # [7] : type (SNV, Deletion, Insertion)
    # [8] : genotype (het., hom. ALT)
    # [9] : base position in sequence file
    # [10] : signal position

    #chromosome number, nucleotide position-1, nucleotide position
    #bed_format = (variant[0], variant[1]-1, variant[1])
    #snp index, chromosome number, nucleotide position, 1, REF/ALT
    #snp_format = (i, variant[0], variant[1], 1, f'{variant[3]}/{variant[4]}')

if __name__ == '__main__': 
    main()