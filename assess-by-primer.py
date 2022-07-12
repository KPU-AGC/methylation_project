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
        help='Path to the directory containing input files'
    )
    parser.add_argument(
        '--output',
        dest='output_path',
        action='store',
        default=None,
        type=pathlib.Path,
        help='Output path'
    )

    args = parser.parse_args()
    if args.output_path is None: 
        args.output_path = args.input_path
    return (args.input_path, args.output_path)

def main(): 
    input_path, output_path = parse_args()
    dict_df = dict()

    for csv_path in input_path.glob('*.csv'): 
        sample_name = csv_path.stem.split('.')[0]
        with open(csv_path, 'r') as sample_file: 
            csv_reader = csv.reader(sample_file)
            dict_df[sample_name] = [x for x in csv_reader]
        #dict_df[sample_name] = pd.read_csv(csv, header=0, index_col=False) 

    sample_names = list(dict_df.keys())
    primer_names = [x[0] for x in dict_df[sample_names[0]]]
    primer_names.remove('primer')
    print(sample_names)
    print(primer_names)

    data_by_primer = dict()

    for primer_name in primer_names: 
        list_data = []
        for sample_name in dict_df.keys(): 
            for line in dict_df[sample_name]: 
                if line[0] == primer_name: 
                    data = line
                    data.insert(1, sample_name)
                    list_data.append(data)
                    break
        data_by_primer[primer_name] = sorted(list_data, key=lambda sample: sample[1])
    
    for primer_name in data_by_primer.keys(): 
        csv_output_path = output_path.joinpath(f'{primer_name}_summary.csv')
        with open(csv_output_path, 'w', newline='') as csv_file: 
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(
                (
                    'primer_name',
                    'sample_name',
                    '#CG_total',
                    '#CG_covered',
                    'mean_methylation',
                    'std_methylation',
                    'mean_coverage',
                    'std_coverage',
                    '#CG<mean',
                )
            )
            csv_writer.writerows(data_by_primer[primer_name])

    '''
    primer_names = list(dict_df[sample_names[0]].loc[:,'primer'])

    data_by_primer = []

    for primer_name in primer_names: 
        list_data = []
        for sample_name in dict_df.keys(): 
            data = dict_df[sample_name].loc[dict_df[sample_name]['primer']==primer_name]
            print(data)
            #data.insert(0, primer_name)
            #data.insert(1, sample_name)
            #list_data.append(data)
        data_by_primer.append(list_data)
    
    print(data_by_primer)
    '''

if __name__ == '__main__': 
    main()