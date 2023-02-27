import argparse
import pathlib
import pandas as pd
import csv

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

    data = dict()

    #Import sequence counts if possible
    if args.nom_path: 
        seq_counts = dict()
        with open(args.nom_path, 'r', encoding='utf-8-sig') as seq_file: 
            csv_reader = csv.reader(seq_file)
            for line in csv_reader: 
                seq_counts[line[0]]=int(line[1].strip('\n'))
        
        print(seq_counts)

    #Populate CSV paths with csvs for each site
    for csv_path in args.input_path.glob('*.csv'): 
        sample_name = csv_path.stem.split('_')[0]
        primer_name = csv_path.stem.split('_')[-1].strip('processed_summary.csv')
        
        if primer_name not in data: 
            data[primer_name] = []

        #Calculations
        csv_data = pd.read_csv(csv_path)
        mean_methylation = csv_data['methylation_percentage'].mean()
        std_methylation = csv_data['methylation_percentage'].std()
        mean_coverage = csv_data['coverage'].mean()
        std_coverage = csv_data['coverage'].std()
        
        row_data = [
            sample_name,
            mean_methylation,
            std_methylation,
            mean_coverage,
            std_coverage,
        ]
        
        if args.nom_path: 
            row_data.append((mean_coverage/seq_counts[sample_name])*30000)
            row_data.append((std_coverage/seq_counts[sample_name])*30000)

        data[primer_name].append(row_data)

    #Output
    for primer in data: 
        output_path = args.output_path.joinpath(f'{primer}_summary.csv')
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

if __name__ == '__main__': 
    main()