'''
process-bismark-coverage.py - script to convert coverage files into organized human-readable .csv files

Functions
---------
import_primer_data : import primer data into named tuple 
get_coverage_paths : if coverage_path argument is a directory, get all coverage file paths
import_coverage_data : import coverage data into DataFrame
process_coverage_data : generate individual and summary data
get_off_target_df : report methylation calls that do not map to target sites
output_summary_data : output summary .csv files
output_primer_df : output individual primer .csv files
'''
__author__ = 'Michael Ke'
__version__ = '1.0.0'
__comments__ = 'Stable.'

#Standard libraries
import argparse
import pathlib
import csv
from collections import namedtuple
#Third-party module
import pandas as pd
import numpy as np

def import_primer_data(primer_path: pathlib.Path) -> list:
    '''Import primer data as a list of namedtuples'''

    Primer = namedtuple(
        'Primer',
        'chromosome, start, end, primer, sequence, num_cg'
    )

    primer_data = []

    with open(primer_path, 'r', encoding='utf-8-sig') as primer_file: 
        csv_reader = csv.reader(primer_file, delimiter='\t')

        for primer in csv_reader: 
            num_cg = primer[4].count('CG')
            primer_data.append(
                    Primer(
                        primer[0],
                        primer[1],
                        primer[2],
                        primer[3],
                        primer[4],
                        num_cg
                )
            )
    
    return primer_data

def get_coverage_paths(coverage_dir_path: pathlib.Path) -> list: 
    '''Create a list of paths for the coverage directory'''
    coverage_paths = []

    for coverage_path in coverage_dir_path.glob('*.cov'): 
        coverage_paths.append(coverage_path)
    
    if not(coverage_paths): 
        print('ERROR: no coverage paths detected')
    else: 
        return coverage_paths

def import_coverage_data(coverage_path: pathlib.Path) -> pd.DataFrame: 
    '''import coverage data into a DataFrame'''
    
    column_names = (
        'chromosome',
        'start',
        'end',
        'methylation_percentage',
        'count_methylated',
        'count_unmethylated',
    )

    coverage_df = pd.read_table(
        coverage_path, 
        names=column_names,
        dtype={'chromosome':'str'}
    )

    return coverage_df

def process_coverage_data(primer_data: list, coverage_df: pd.DataFrame): 
    #Summary data
    summary_data = []
    primer_dfs = dict()

    #Perform for each primer set
    for primer in primer_data:

        #Get Dataframe of CG positions within primer region 
        region_df = coverage_df[
            (coverage_df.chromosome == primer.chromosome) 
            & (coverage_df.start > int(primer.start)) 
            & (coverage_df.end < int(primer.end))
        ]

        if region_df.empty:
            #print(f'No data for {primer.primer}')
            pass
        else:
            #For each CG position, get the position relative to the sequence file 
            # in 0-based coordinate, then store in dataframe
            region_df.loc[:,'seq_pos'] = region_df.loc[:,'start'] - 1 - int(primer.start)

            #Sum unmethylated and methylated read calls to get total coverage
            region_df.loc[:,'coverage'] = (
                region_df.loc[:,'count_methylated'] 
                + region_df.loc[:,'count_unmethylated']
            )

            #TODO: this function throws an error that's pretty aids. Need to fix in the future.
            #Determine basecall for each CG position (is it C or G)
            basecalls = []
            for cg in region_df.iterrows(): 
                basecalls.append(primer.sequence[cg[1]['seq_pos']])
            region_df.loc[:,'basecall'] = basecalls

            #Check if the current primer set is on the positive or negative strand
            strand = True
            if (
                (primer.primer.split('-')[1][0:3] == 'BSN') 
                or (primer.primer.split('-')[1][0:2] == 'BN')
            ):  
                strand = False
            
            #Generate primer region DataFrame depending on assigned basecall
            if strand is True: 
                primer_df = region_df.loc[region_df['basecall']=='C']
            else: 
                primer_df= region_df.loc[region_df['basecall']=='G']

            #Summary statistics
            number_cg = len(primer_df)
            methylation_avg = primer_df['methylation_percentage'].mean()
            methylation_std = primer_df['methylation_percentage'].std()
            if np.isnan(methylation_std): 
                methylation_std = 0
            coverage_avg = primer_df['coverage'].mean()
            coverage_std = primer_df['coverage'].std()
            if np.isnan(coverage_std):
                coverage_std = 0
            under_coverage = sum(primer_df['coverage'] < coverage_avg)

            summary_data.append(
                (
                    primer.primer,
                    primer.num_cg,
                    number_cg,
                    methylation_avg, 
                    methylation_std, 
                    coverage_avg, 
                    coverage_std,
                    under_coverage
                    )
                )

            primer_dfs[primer.primer] = primer_df

            #print(primer_df)
            #print(f'{primer.primer} mean methylation: {methylation_avg}')
            #print(f'{primer.primer} std methylation: {methylation_std}')
            #print(f'{primer.primer} mean coverage: {coverage_avg}')
            #print(f'{primer.primer} std coverage: {coverage_std}')
            #print(f'{primer.primer} number CG under mean coverage: {under_coverage}')

    return summary_data, primer_dfs

def get_off_target_df(primer_data: list, coverage_df: pd.DataFrame) -> pd.DataFrame: 
    '''Figure out what other sequences remain.'''
    off_target_df = coverage_df
    #print(off_target_df)
    for primer in primer_data: 
        off_target_df.drop(
            off_target_df[
                (off_target_df.chromosome == primer.chromosome)
                & (off_target_df.start > int(primer.start))
                & (off_target_df.start < int(primer.end))
            ].index,
            inplace=True,
        )

    return off_target_df

def output_summary_data(
        output_path: pathlib.Path, 
        summary_data: list, 
        sample_name: str,
    ) -> None: 
    '''Output summary data for coverage file to .csv'''
    
    summary_names=(
        'primer',
        '#CG_total',
        '#CG_covered',
        'mean_methylation',
        'std_methylation',
        'mean_coverage',
        'std_coverage',
        '#CG<mean',
    )

    summary_df = pd.DataFrame(summary_data, columns=summary_names)

    output_csv_path = output_path.joinpath(f'{sample_name}_summary.csv')

    summary_df.to_csv(output_csv_path, index=False)

def output_primer_df(
        output_path: pathlib.Path, 
        primer_df: pd.DataFrame, 
        primer: str,
        sample_name: str,
    ) -> None:
    '''Output primer data to a csv.'''
    primer_data_path = output_path.joinpath(f'{sample_name}_{primer}_processed.csv')
    primer_df.to_csv(primer_data_path, index=False)

def parse_args(): 
    parser = argparse.ArgumentParser("Script to analyze the coverage files")
    parser.add_argument(
        'coverage_path',
        action='store', 
        type=pathlib.Path,
        help='Path to coverage file(s) generated by bismark methylation extractor.'
    )
    parser.add_argument(
        'primer_path',
        action='store',
        type=pathlib.Path,
        help='Path to primer file.'
    )
    parser.add_argument(
        '-o',
        dest='output_path',
        action='store',
        default=None,
        type=pathlib.Path,
        help='Directory to output to.'
    )
    args = parser.parse_args()
    if args.output_path is None: 
        args.output_path = args.coverage_path.parent
    return args

def main(): 
    args = parse_args()

    #NOTE:
    #THE PRIMER FILE HAS COORDINATES IN 0-BASED, HALF-OPEN. 
    #THE COVERAGE REPORT IS IN 1-BASED, FULLY-CLOSED.
    #THEREFORE, TO MAINTAIN CONSISTENCY IN MATH, WE WILL CONVERT THE COVERAGE REPORT
    #COORDINATES TO 0-BASED AS WELL. 
    primer_data = import_primer_data(args.primer_path)
    
    #MUTE chained assignment warning
    pd.options.mode.chained_assignment = None

    if args.coverage_path.is_dir(): 
        coverage_paths = get_coverage_paths(args.coverage_path)

        summary_dir = args.output_path.joinpath('summaries')
        primer_dir = args.output_path.joinpath('individual_primer')
        summary_dir.mkdir(exist_ok=True)
        primer_dir.mkdir(exist_ok=True)
        print(summary_dir)
        print(primer_dir)

        for coverage_path in coverage_paths: 
            sample_name = coverage_path.stem.split('_')[0]
            coverage_df = import_coverage_data(coverage_path)
            summary_data, primer_dfs = process_coverage_data(primer_data, coverage_df)
            off_target_df = get_off_target_df(primer_data, coverage_df)

            for primer in primer_dfs: 
                output_primer_df(
                    primer_dir, 
                    primer_dfs[primer], 
                    primer,
                    sample_name,
                )

            output_summary_data(
                summary_dir,
                summary_data,
                sample_name
            )

            #Output remainder mapped
            output_primer_df(
                primer_dir,
                off_target_df,
                'off-target',
                sample_name,
            )

    elif args.coverage_path.is_file(): 
        #Example args.coverage_path name: 
        #BSX-18-B-None-Multiplex_S4_L001_1_bismark_bt2_pe.nonCG_filtered.bismark.cov
        sample_name = args.coverage_path.stem.split('_')[0]
        coverage_df = import_coverage_data(args.coverage_path)
        summary_data, primer_dfs = process_coverage_data(primer_data, coverage_df)
        off_target_df = get_off_target_df(primer_data, coverage_df)

        for primer in primer_dfs: 
            output_primer_df(
                args.output_path,
                primer_dfs[primer], 
                primer,
                sample_name,
            )

        output_summary_data(
            args.output_path,
            summary_data,
            sample_name
        )

        #Output remainder mapped
        output_primer_df(
            args.output_path,
            off_target_df,
            'off-target',
            sample_name,
        )

    else:
        print('Path is neither a directory or file.')
        exit()

if __name__ == '__main__': 
    main()