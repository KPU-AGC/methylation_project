"""
cgvisualizer.py - generate % methylation plots for each sample

General workflow
----------------
1) Determine if coverage_file_path is a directory or path; loop through all coverage_files
   if it is a directory
2) Import primer data into dict
3) Import coverage data for sample, and separate by primer set (target sites) and record
   data + summary statistics into dictionary
4) Create a figure for the sample where it plots: 
    a) % methylation at each CG position
    b) % of coverage against max coverage value at each CG position
   Plot is generated for each target site, and combined into one figure. 
   Target sites with no data are marked. 

Functions
---------
import_coverage_data : from coverage file, import all of the data into a DataFrame
import_primer_data : from primer.bed file, import all primer data into a dictionary
get_region_data : split coverage DataFrame into list of dictionaries, where each dictionary
                  contains methylation information about that region
generate_figure : generate figure containing methylation plots for each target site (region)
"""

__author__ = 'Michael Ke'
__version__ = '1.0.0'
__comments__ = 'Should be stable.'

#Standard libraries
from statistics import mean
import argparse
import csv
import pathlib
import re
from math import ceil
#Third-party modules
from matplotlib import pyplot
import numpy as np
import pandas as pd

def import_coverage_data(path: pathlib.Path) -> pd.DataFrame:
    '''Return pandas table with coverage file data'''
    return pd.read_table(
        path, 
        names=(
            'chromosome',
            'start',
            'end',
            'methylation_percentage',
            'count_methylated',
            'count_unmethylated',
            ), 
        dtype={'chromosome':'str'},
        )

def import_primer_data(path: pathlib.Path) -> dict: 
    """Return primer data from bed file specifiying target regions"""

    primer_data = dict()

    with open(path, 'r', encoding='utf-8-sig') as primer_file: 
        csv_reader = csv.reader(
            primer_file, 
            delimiter='\t'
        )

        for primer in csv_reader: 
            num_cg = primer[4].count('CG')

            cg_iter = re.finditer('CG', primer[4])
            cg_locations = []

            for cg in cg_iter:
                #1-based index 
                cg_locations.append(cg.start()+int(primer[1])+1)
            
            primer_data[primer[3]] = {
                'chromosome':primer[0],
                'start':primer[1],
                'end':primer[2],
                'sequence':primer[4],
                'num_cg':num_cg, 
                'cg_pos':cg_locations,
            }

    return primer_data   

def get_region_data(primer_data: dict, coverage_data: pd.DataFrame) -> list:
    """
    Return a list of regions and their methylation calls at each CG site

    Parameters
    ----------
    primer_data : list
        contains information about each target site
    coverage_data : pd.DataFrame
        all methylation calls for the sample

    Return
    ------
    List of regions and their methylation calls at each CG site
    """

    data = []

    for primer in primer_data: 
        chromosome=primer_data[primer]['chromosome']
        start=primer_data[primer]['start']
        end=primer_data[primer]['end']        
        sequence=primer_data[primer]['sequence']

        #Get Dataframe of CG positions within primer region 
        region_df = coverage_data[
            (coverage_data.chromosome == chromosome) 
            & (coverage_data.start > int(start)) 
            & (coverage_data.end < int(end))
            ]
        
        #For each CG position, get the position relative to the sequence file 
        # in 0-based coordinate, then store in dataframe
        region_df.loc[:,'seq_pos'] = region_df.loc[:,'start'] - 1 - int(start)

        #Determine basecall for each CG position (is it C or G)
        basecalls = []
        for cg in region_df.iterrows(): 
            basecalls.append(sequence[cg[1]['seq_pos']])
        region_df.loc[:,'basecall'] = basecalls

        #Check if the current primer set is on the positive or negative strand
        strand = True
        if (
            (primer.split('-')[1][0:3] == 'BSN') 
            or (primer.split('-')[1][0:2] == 'BN')
        ):  
            strand = False

        #Generate DataFrame with all of the data for that specific primer set   
        if strand is True: 
            primer_df = region_df.loc[region_df['basecall']=='C']
        else: 
            primer_df= region_df.loc[region_df['basecall']=='G']

        primer_df.loc[:,'coverage'] = (
            primer_df.loc[:, 'count_methylated'] 
            + primer_df.loc[:,'count_unmethylated']
        )

        region_data = {
            'primer':primer,
            'positions':primer_df.loc[:,'start'].to_list(),
            'methylation':primer_df.loc[:,'methylation_percentage'].to_list(),
            'coverage':primer_df.loc[:,'coverage'].to_list(),
        }
        try: 
            max_cov = max(region_data['coverage'])
            region_data['relative_coverage'] = [(cov/max_cov)*100 for cov in region_data['coverage']]
            print(region_data['relative_coverage'])
        except ValueError: 
            print('No data')

        data.append(region_data)

    return data

def generate_figure(
        primer_data: dict, 
        region_data: list, 
        sample_name: str,
        output_path: pathlib.Path,
        show_fig_flag: bool,
    ) -> None:
    """
    Generate a figure containing methylation plots for each target site. 

    Parameters
    ----------
    primer_data : dict
        key (primer set) yields information about the primer set
    region_data : list
        each entry corresponds to methylation data for target site
    sample_name : str
        name of sample that region_data belongs to
    output_path : pathlib.Path
        directory to output figures (.png) to
    show_fig_flag : bool
        if set to True, will preview the figure
    
    Returns
    -------
    None; writes figure .png to specified directory. 
    
    """
    def _get_fig_size(num_primers: int) -> int: 
        fig_height = ceil(num_primers/3)
        return fig_height, 3
    
    #Calculate the number of subplots that should exist
    #Currently, our primer sets look like this: 
    #BS-Run-1 : 30 sets
    #BS-Run-2 : 20 sets
    #BSX-Run-2: 16 sets
    #BSX-Run-3: 21 sets

    #TODO: do not hardcode these plots
    #Current implementation - based on the number of primers, generate 
    # the appropriate width and height params to fit all samples
    num_primers = len(primer_data)
    #print(num_primers)
    #if num_primers == 30: 
    #    height = 10
    #    width = 3
    #elif num_primers == 21:
    #    height = 7
    #    width = 3
    #elif num_primers == 20: 
    #    height = 10
    #    width = 2
    #elif num_primers == 17: 
    #    height = 6
    #    width = 3
    #else: 
    #    print('NUMBER OF PRIMERS NOT USED IN BOVITEQ PROJECT.')
    #    print('PROGRAM CANNOT CONTINUE. FUTURE IMPLEMENTATION WILL FIX THIS ISSUE')
    #    exit() 

    height, width = _get_fig_size(num_primers)

    fig, axs = pyplot.subplots(height, width, figsize=(8.5,11))
    
    region_index = 0
    #Iterate through each of the 21 plots
    #Each plot specified as axs[i][j]
    for i in range(height): 
        for j in range(width): 
            if region_index < num_primers:
                primer = region_data[region_index]['primer']
                #pos_data = region_data[region_index]['positions']
                methyl_data = region_data[region_index]['methylation']
                #cg_pos = primer_data[primer]['cg_pos']

                if methyl_data: 
                    #Assign colours - mean methylation
                    #Green = Methylated, >70%
                    #Yellow = 30-70%
                    #Blue = Unmethylated, <30%
                    
                    max_coverage = max(region_data[region_index]['coverage'])

                    mean_methylation = mean(methyl_data)
                    if mean_methylation >= 70.0: 
                        colour = 'forestgreen'
                    elif  30 < mean_methylation < 70.0: 
                        colour = 'wheat'
                    else: 
                        colour = 'lightskyblue'

                    axs[i][j].bar(
                        range(len(methyl_data)),
                        methyl_data, 
                        width=1,
                        align='center',
                        color=colour,
                    )
                    axs[i][j].plot(
                        region_data[region_index]['relative_coverage'],
                        'red'
                        )
                    axs[i][j].text(
                        0,
                        3,
                        f'max(coverage)={str(max_coverage)}',
                        fontsize='x-small',
                    )
                
                else: 
                    axs[i][j].bar(
                        0,
                        0,
                        width=1,
                        align='edge',
                    )
                    axs[i][j].text(
                        0.5,
                        50,
                        'No data',
                        horizontalalignment='center',
                        verticalalignment='center',  
                    )
                axs[i][j].set(
                    ylim=[0,100],
                    title=f'{primer}',
                )
                axs[i][j].set_yticks(
                    (0, 50.0, 100.0),
                )
                axs[i][j].set_yticklabels(
                    ('0%', '50%', '100%')
                )
            else: 
                pass
            region_index = region_index + 1
    
    pyplot.tight_layout()
    fig_output_path = output_path.joinpath(f'{sample_name}.png')
    fig.subplots_adjust(top=0.94)
    fig.suptitle(sample_name)
    #pyplot.tight_layout()

    fig.savefig(fig_output_path, transparent=False, dpi=300)

    if show_fig_flag:
        pyplot.show()

def parse_args(): 
    parser = argparse.ArgumentParser("Program to generate images of CG island")

    parser.add_argument(
        'coverage_file_path',
        type=pathlib.Path,
        action='store',
        help='Path to coverage file(s)'
    )
    parser.add_argument(
        'primer_file_path',
        type=pathlib.Path,
        action='store',
        help='Path to .csv file with primer information.',
    )
    parser.add_argument(
        '-o',
        '--output',
        dest='output_path',
        type=pathlib.Path,
        action='store',
        default=None,
        help='Directory to output results to.'
    )
    parser.add_argument(
        '-v',
        '--verbose',
        dest='show_fig_flag',
        action='store_true',
        help='Show each figure as they are being generated'
    )
    parser.add_argument(
        '-s',
        '--snpsplit',
        dest='snp_flag',
        action='store_true',
        help='Set if you are processing SNPsplit BAMs.',
    )

    args = parser.parse_args()
    
    if not args.output_path: 
        args.output_path = args.coverage_file_path

    return args

def main(): 
    args = parse_args()

    primer_data = import_primer_data(args.primer_file_path)

    if args.coverage_file_path.is_dir(): 
        for coverage_file_path in args.coverage_file_path.glob('*.cov'): 
            #Define correct sample ID name
            if args.snp_flag: 
                sample_name = coverage_file_path.stem.split('_')[0]
                pool_name = coverage_file_path.stem.split('.')[-2]
                sample_id = f'{sample_name}_{pool_name}'
            else: 
                sample_id = coverage_file_path.stem.split('_')[0]

            coverage_data = import_coverage_data(coverage_file_path)
            region_data = get_region_data(primer_data, coverage_data)
            generate_figure(
                primer_data, 
                region_data, 
                sample_id, 
                args.output_path,
                args.show_fig_flag,
            )
    elif args.coverage_file_path.is_file(): 
        sample_name = args.coverage_file_path.stem.split('_')[0]
        coverage_data = import_coverage_data(args.coverage_file_path)
        region_data = get_region_data(primer_data, coverage_data)
        generate_figure(
            primer_data, 
            region_data, 
            sample_name, 
            args.output_path,
            args.show_fig_flag,
        )
    else: 
        print('Coverage path is neither a directory or file.')
        exit()

if __name__ == '__main__': 
    main()