"""
assess-by-primer.py - script for calculating coverage box plots for each primer site. 

Functions
---------
import_seq_counts : from .csv file, get number of post-processing sequences for each sample
import_bisulfite_data : import data from sample-primer.csv output of process-ngs-coverage.
output_primer_summaries : output .csv files for each primer with summary statistics
generate_figure : generate box plot image

"""
from matplotlib import pyplot
import argparse
import pathlib
import pandas as pd
import numpy as np
import csv

def import_seq_counts(nom_path: pathlib.Path) -> dict: 
    '''Import the total sequence counts after mapping and filtering per sample.'''
    seq_counts = dict()

    with open(nom_path, 'r', encoding='utf-8-sig') as seq_file: 
        csv_reader = csv.reader(seq_file)
        for line in csv_reader: 
            seq_counts[line[0]]=int(line[1].strip('\n'))
    
    return seq_counts

def import_bisulfite_data(input_dir: pathlib.Path, seq_counts: dict, num_nom: int) -> dict: 

    data = dict()
    data_df = dict()

    labels = [
        'sample_name',
        'mean_methylation',
        'std_methylation',
        'mean_coverage',
        'std_coverage',
        'norm_coverage',
        'norm_std_coverage',
    ]  

    #Initial data import
    for csv_path in input_dir.glob('*.csv'): 

        sample_name = csv_path.stem.split('_')[0]
        primer_name = csv_path.stem.split('_')[-1].strip('processed_summary.csv')
        
        #Create new list with primer name if not present
        if primer_name not in data: 
            data[primer_name] = []

        #Calculations
        csv_data = pd.read_csv(csv_path)
        mean_methylation = csv_data['methylation_percentage'].mean()
        std_methylation = csv_data['methylation_percentage'].std()
        mean_coverage = csv_data['coverage'].mean()
        std_coverage = csv_data['coverage'].std()
        if seq_counts:
            norm_coverage = (mean_coverage/seq_counts[sample_name])*num_nom
            norm_std_coverage = (std_coverage/seq_counts[sample_name])*num_nom
        else: 
            norm_coverage = mean_coverage/1
            norm_std_coverage = std_coverage/1

        #Add data to appropriate list in dictionary (based on primer)
        row_data = [
            sample_name,
            mean_methylation,
            std_methylation,
            mean_coverage,
            std_coverage,
            norm_coverage,
            norm_std_coverage,
        ]

        data[primer_name].append(row_data)

    #Conversion to DataFrames
    for primer in data: 
        data_df[primer] = pd.DataFrame(data[primer], columns=labels)

    return data_df

def output_primer_summaries(data: dict, output_dir: pathlib.Path, tag: str) -> None: 
    '''Output primer summaries for each primer DataFrame'''
    
    #Per-primer summary
    for primer in data: 
        #Generate output path
        if tag: 
            output_path = output_dir.joinpath(f'{tag}_{primer}_summary.csv')
        else: 
            output_path = output_dir.joinpath(f'{primer}_summary.csv')
        #Output
        with open(output_path, 'w', newline='') as output_file: 
            data[primer].to_csv(
                output_file, 
                header=True,
                index=False,
            )

def generate_figure(data: dict, output_path: pathlib.Path, tag: str) -> None: 
    '''Create box plot'''
    
    fig, axe = pyplot.subplots()
    
    id_list = []
    data_list = []

    #Add data to a list
    for primer in data: 
        norm_coverage_series = data[primer]['norm_coverage']
        id_list.append(f'{primer} (n={str(norm_coverage_series.size)})')
        data_list.append(norm_coverage_series)

    #Layout parameters
    axe.boxplot(data_list, labels=id_list)
    pyplot.xticks(rotation='vertical')
    pyplot.tight_layout()

    #pyplot.show()

    #Generate output path
    if tag:
        fig_path = output_path.joinpath(f'{tag}_primer_boxplot_summary.png')
    else:
        fig_path = output_path.joinpath(f'primer_boxplot_summary.png')

    fig.savefig(fig_path, transparent=False, dpi=300)

def parse_args(): 
    parser = argparse.ArgumentParser("Assess summary statistics by primer")
    parser.add_argument(
        "input_path",
        action='store', 
        type=pathlib.Path,
        help='Path to the directory containing .csv files from process-bismark-coverage.py'
    )
    parser.add_argument(
        '-o', 
        '--output',
        dest='output_path',
        action='store',
        default=None,
        type=pathlib.Path,
        help='Output path'
    )
    parser.add_argument(
        '-n', 
        '--nom_file',
        dest='nom_path',
        action='store',
        default=None,
        type=pathlib.Path,
        help='Path to file containing sample names and final mapped sequence count after filtering.'
    )
    parser.add_argument(
        '-t',
        '--tag',
        dest='tag',
        default='',
        type=str,
        help='Optional prefix to all files generated.'
    )
    parser.add_argument(
        '-c',
        '--nom_count',
        dest='num_nom',
        action='store',
        default=30000,
        type=int,
        help='Number of reads to normalize to. Default=30000'
    )

    args = parser.parse_args()

    if not(args.output_path): 
        args.output_path = args.input_path
    
    return args

def main(): 
    args = parse_args()

    #Import sequence counts if possible
    if args.nom_path: 
        seq_counts = import_seq_counts(args.nom_path)
    else: 
        seq_counts = []

    #Populate CSV paths with csvs for each site
    data = import_bisulfite_data(args.input_path, seq_counts, args.num_nom)

    #Output - DataFrame outputs
    output_primer_summaries(data, args.output_path, args.tag)

    generate_figure(data, args.output_path, args.tag)

if __name__ == '__main__': 
    main()