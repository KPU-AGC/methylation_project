from matplotlib import pyplot
import argparse
import pathlib
import pandas as pd
import numpy as np
import csv

def import_seq_counts(nom_path: pathlib.Path): 
    '''Import the total sequence counts after mapping and filtering per sample.'''
    seq_counts = dict()

    with open(nom_path, 'r', encoding='utf-8-sig') as seq_file: 
        csv_reader = csv.reader(seq_file)
        for line in csv_reader: 
            seq_counts[line[0]]=int(line[1].strip('\n'))
    
    return seq_counts

def import_bisulfite_data(input_dir: pathlib.Path, seq_counts): 

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
            norm_coverage = (mean_coverage/seq_counts[sample_name])*30000
            norm_std_coverage = (std_coverage/seq_counts[sample_name])*30000

        #Add data to appropriate list in dictionary (based on primer)
        row_data = [
            sample_name,
            mean_methylation,
            std_methylation,
            mean_coverage,
            std_coverage,
        ]

        if seq_counts: 
            row_data.append(norm_coverage)
            row_data.append(norm_std_coverage)

        data[primer_name].append(row_data)

    #Conversion to DataFrames
    for key in data: 
        data_df[key] = pd.DataFrame(data[key], columns=labels)

    return data_df

def output_primer_summaries(data: dict, output_dir: pathlib.Path): 
    '''Output primer summaries for each primer DataFrame'''
    
    for primer in data: 
        output_path = output_dir.joinpath(f'{primer}_summary.csv')
        with open(output_path, 'w', newline='') as output_file: 
            data[primer].to_csv(
                output_file, 
                header=True,
                index=False,
            )

def output_primer_summaries_old(data: dict, output_path: pathlib.Path): 
    '''Output primer summaries'''
    for primer in data: 
        output_path = output_path.joinpath(f'{primer}_summary.csv')
        with open(output_path, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(
                (
                    'sample_name',
                    'mean_methylation',
                    'std_methylation',
                    'mean_coverage',
                    'std_coverage',
                    'normalized_mean_coverage',
                    'normalized_std_coverage',
                )
            )
            csv_writer.writerows(data[primer])

def generate_figure(data: dict, output_path: pathlib.Path): 
    ''''''
    fig, axe = pyplot.subplots()
    
    id_list = []
    data_list = []

    #Add data to a list
    for key in data: 
        norm_coverage_series = data[key]['norm_coverage']
        id_list.append(f'{key} (n={str(norm_coverage_series.size)})')
        data_list.append(norm_coverage_series)

    axe.boxplot(data_list, labels=id_list)
    pyplot.xticks(rotation='vertical')
    pyplot.tight_layout()
    #pyplot.show()
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
    data = import_bisulfite_data(args.input_path, seq_counts)

    #Output - DataFrame outputs
    output_primer_summaries(data, args.output_path)

    generate_figure(data, args.output_path)

if __name__ == '__main__': 
    main()