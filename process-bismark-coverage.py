import argparse
import pathlib
import pandas as pd
import csv
from collections import namedtuple

class CoverageHandler(): 
    pass

def parse_args(): 
    parser = argparse.ArgumentParser("Script to analyze the coverage files")
    parser.add_argument(
        'coverage_path',
        action='store', 
        type=pathlib.Path,
        help='Path to coverage files'
    )
    parser.add_argument(
        'primer_path',
        action='store',
        type=pathlib.Path,
        help='Path to primer file'
    )
    parser.add_argument(
        '-o',
        dest='output_path',
        action='store',
        default=None,
        type=pathlib.Path,
        help='Output path'
    )
    args = parser.parse_args()
    if args.output_path is None: 
        args.output_path = args.coverage_path.parent
    return (args.coverage_path, args.primer_path, args.output_path)

def main(): 
    coverage_path, primer_path, output_path = parse_args()
    coverage_df = pd.read_table(coverage_path, names=('chromosome','start','end','methylation_percentage','count_methylated','count_unmethylated'))
    #Import primers
    Primer = namedtuple('Primer','chromosome, start, end, primer, sequence')
    with open(primer_path, 'r', encoding='utf-8-sig') as primer_file: 
        csv_reader = csv.reader(primer_file, delimiter='\t')
        primer_data = []
        for primer in map(Primer._make, csv_reader): 
            primer_data.append(primer)

    #Perform for each primer set
    for primer in primer_data: 
        primer_df = coverage_df[(coverage_df.chromosome == primer.chromosome) & (coverage_df.start > int(primer.start)) & (coverage_df.end < int(primer.end))]
        #print(primer_df)
        print(primer.sequence)

if __name__ == '__main__': 
    main()