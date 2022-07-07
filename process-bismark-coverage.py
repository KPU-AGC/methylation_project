import argparse
import pathlib
import pandas as pd
import csv

class CoverageHandler(): 
    pass

def parse_args(): 
    parser = argparse.ArgumentParser("Script to analyze the coverage files")
    parser.add_argument(
        'coverage_file',
        action='store', 
        type=pathlib.Path,
        help='Path to coverage files'
    )
    parser.add_argument(
        'primer_file',
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
        args.output_path = args.coverage_file.parent
    return (args.coverage_file, args.primer_file, args.output_path)

def main(): 
    coverage_path, primer_path, output_path = parse_args()

if __name__ == '__main__': 
    main()