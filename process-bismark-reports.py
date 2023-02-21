#standard libraries
import argparse
import pathlib
import csv
import logging

def parse_args(): 
    parser = argparse.ArgumentParser('Program to combine all of the bismark reports')
    parser.add_argument(
        'mapped_bam_path',
        type=pathlib.Path,
        action='store',
        help='Path to directory containing all reports.'
    )
    parser.add_argument(
        '-o', '--output',
        type=pathlib.Path,
        dest='output_path',
        default=None,
        action='store',
        help='Output directory.'
    )
    args = parser.parse_args()
    if not args.output_path: 
        args.output_path = args.input_path
    
    return args

def main(): 
    args = parse_args()

    report_paths = []

    for file_path in args.mapped_bam_path.glob('*.txt'): 
        report_paths.append(file_path)

    #Output all paths of files that are being processed
    print(report_paths)

    qc_data = []

    for report_path in report_paths: 
        with open(report_path, 'r') as report_file: 
            lines = report_file.readlines()
            total_seq = lines[6].split('\t')[1].strip('\n')
            mapped = lines[7].split('\t')[1].strip('\n')
            efficiency = lines[8].split('\t')[1].strip('\n')
            unmapped = lines[9].split('\t')[1].strip('\n')
            sample_name = report_path.stem.split('_')[0]
            qc_data.append(
                (
                    sample_name, 
                    total_seq, 
                    mapped, 
                    unmapped, 
                    efficiency
                )
            )
    
    output_path = args.output_path.joinpath('summary_bismark.csv')
    with open(output_path, 'w', newline='') as csv_file: 
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(
            (
                'sample_name',
                'total_seq',
                'mapped',
                'unmapped',
                'efficiency',
            )
        )
        csv_writer.writerows(qc_data)

if __name__ == '__main__': 
    main()