"""

"""
from matplotlib import pyplot
import numpy as np
from statistics import mean
import argparse
import pandas as pd
import csv
import pathlib
import re
from collections import namedtuple

def import_coverage_data(path: pathlib.Path):
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

def import_primer_data_old(path: pathlib.Path): 
    '''Return primer data from bed file specifiying target regions'''
    Primer = namedtuple(
        'Primer',
        'chromosome, start, end, primer, sequence, num_cg, cg_locations',
        )
    primer_data = []

    with open(path, 'r', encoding='utf-8-sig') as primer_file: 
        csv_reader = csv.reader(primer_file, delimiter='\t')
        primer_data = []

        for primer in csv_reader: 
            num_cg = primer[4].count('CG')

            cg_iter = re.finditer('CG', primer[4])
            cg_locations = []
            for cg in cg_iter: 
                cg_locations.append(cg.start())

            primer_data.append(Primer(
                primer[0],
                primer[1],
                primer[2],
                primer[3],
                primer[4],
                num_cg,
                cg_locations
            ))

    return primer_data

def import_primer_data(path: pathlib.Path): 
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

def get_region_data(primer_data: dict, coverage_data: pd.DataFrame):
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


def generate_figure(primer_data: dict, region_data: list, coverage_path: pathlib.Path):

    fig, axs = pyplot.subplots(7, 3, figsize=(8.5,11))
    
    region_index = 0
    #Iterate through each of the 21 plots
    for i in range(7): 
        for j in range(3): 
            
            primer = region_data[region_index]['primer']
            pos_data = region_data[region_index]['positions']
            methyl_data = region_data[region_index]['methylation']
            cg_pos = primer_data[primer]['cg_pos']

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
                    align='edge',
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
            region_index = region_index + 1
    pyplot.tight_layout()

    sample_name = coverage_path.stem.split('_')[0]
    output_path = coverage_path.parent.joinpath(f'{sample_name}.png')
    fig.subplots_adjust(top=0.94)
    fig.suptitle(sample_name)
    fig.savefig(output_path, transparent=False, dpi=80)
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

    args = parser.parse_args()

    return args

def main(): 
    args = parse_args()

    primer_data = import_primer_data(args.primer_file_path)
    coverage_data = import_coverage_data(args.coverage_file_path)

    region_data = get_region_data(primer_data, coverage_data)
    generate_figure(primer_data, region_data, args.coverage_file_path)

if __name__ == '__main__': 
    main()