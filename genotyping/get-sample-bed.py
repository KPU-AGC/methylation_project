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
        help='Path to the json output of combine-variant-calls.py'
    )
    parser.add_argument(
        'target_sample',
        action='store',
        type=str,
        help='Key of target sample'
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
    return args.json_path, args.target_sample, args.output_path, args.snp_flag, args.het_flag

def variant_output(data, target, output_path): 
    #Create and setup output CSV file
    csv_output_path = output_path.joinpath(f'{target}_variants.bed')
    csv_file = open(csv_output_path, 'w', newline='')
    csv_writer = csv.writer(csv_file, delimiter='\t')
    csv_writer.writerows(data)
    csv_file.close()

def snp_output(data, target, output_path): 
    #Create and setup output CSV file
    csv_output_path = output_path.joinpath(f'{target}_variants.txt')
    csv_file = open(csv_output_path, 'w', newline='')
    csv_writer = csv.writer(csv_file, delimiter='\t')
    csv_writer.writerow(("ID", "Chr", "Position", "SNP_value", "Ref/SNP"))
    csv_writer.writerows(data)
    csv_file.close()


def main(): 
    json_path, target, output_path, snp_flag, het_flag = parse_args()

    json_file = open(json_path, 'r')
    json_data = json.load(json_file)
    json_file.close()
    sample_data = json_data[target]
    bed_variants = set()
    variants = set()
    i = 1
    for variant in sample_data['variants']: 
        variant_include = True
        if snp_flag and not(variant[7] == 'SNV'): 
            variant_include = False
        elif het_flag and not(variant[8] == 'hom. ALT'): 
            variant_include = False
        else: 
            pass
        if variant_include:
            bed_format = (variant[0], variant[1]-1, variant[1])
            snp_format = (i, variant[0], variant[1], 1, f'{variant[3]}/{variant[4]}')
            bed_variants.add(bed_format)
            variants.add(snp_format)
            i = i + 1
    list_variants = sorted(bed_variants)
    list_snps = sorted(variants)
    variant_output(list_variants, target, output_path)
    snp_output(list_snps, target, output_path)

if __name__ == '__main__': 
    main()