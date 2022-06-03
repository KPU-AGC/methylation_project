import argparse
import csv
import pathlib
import json

def parse_args(): 
    parser = argparse.ArgumentParser("Program to combine all of the variant calls into one file")
    parser.add_argument(
        "json_directory",
        action='store', 
        type=pathlib.Path,
        help='Path to the directory containing .json files'
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
        '--no_filter',
        dest='filter_flag',
        action='store_true',
        help='Set to true if you want all variants'
    )
    args = parser.parse_args()
    if args.output_path is None: 
        args.output_path = args.json_directory
    return (args.json_directory, args.output_path, args.filter_flag)

def main(): 
    json_directory, output_path, filter_flag = parse_args()
    #Create and setup output CSV file
    csv_output_path = output_path.joinpath('variants.csv')
    csv_file = open(csv_output_path, 'w', newline='')
    csv_writer = csv.writer(csv_file, delimiter=',')
    csv_writer.writerow(
        (
            'chr',
            'pos',
            'id',
            'ref',
            'alt',
            'qual',
            'filter',
            'type',
            'genotype',
            'basepos',
            'signalpos',
        )
    )

    for json_path in json_directory.glob('*.json'):
        json_file = open(json_path, 'r')
        sample_data = json.load(json_file)
        #HARDCODED: CURRENTLY SAMPLE_GENE-PRIMER-DIRECTION ENFORCED IN JSON NAME
        sample_name = json_path.stem.split('_')[0]
        primer_name = json_path.stem.split('_')[1]
        for variant in sample_data['variants']['rows']:
            #[6] - FILTER
            #[2] - ID - populate with gene-primer
            if (
                variant[6] in ('PASS','MANUAL')
                or filter_flag is True
            ):
                variant[2] = sample_name + '_' + primer_name
                csv_writer.writerow(variant)
        json_file.close()
    csv_file.close()

if __name__ == '__main__':
    main()