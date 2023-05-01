"""
assess-read-depth.py - Determine the depth of reads at each site and then generate plots for each one.

"""

#Standard libraries
import argparse
import pathlib
import csv
from math import ceil

#Third-party modules
from matplotlib import pyplot
import numpy as np
import pandas as pd

def import_primer_data(path: pathlib.Path) -> dict: 
    """Return primer data from bed file specifiying target regions"""

    primer_data = dict()

    with open(path, 'r', encoding='utf-8-sig') as primer_file: 
        csv_reader = csv.reader(
            primer_file, 
            delimiter='\t'
        )

        for primer in csv_reader:
            primer_data[primer[3]] = {
                'chromosome':primer[0],
                'start':primer[1],
                'end':primer[2],
                'sequence':primer[4],
            }

    return primer_data

def import_depth_data(path: pathlib.Path) -> pd.DataFrame: 
    return pd.read_table(
        path,
        sep='\t',
        names=(
            'chromosome',
            'position',
            'count',
        ),
        dtype={'chromsoome':'str'}
    )

def get_region_data(primer_data: dict, depth_data: pd.DataFrame) -> list: 
    """
    Return list of regions with their count data at each base position.

    Parameters
    ----------
    primer_data : list
        contains information about each target site
    depth_data : pd.DataFrame
        all read counts at each base position

    Return
    ------
    List of regions and their methylation calls at each CG site
    """

    data = []

    for primer in primer_data: 
        chromosome=primer_data[primer]['chromosome']
        start=primer_data[primer]['start']
        end=primer_data[primer]['end']     
    
        region_df = depth_data[
            (depth_data.chromosome == chromosome) 
            & (depth_data.position > int(start)) 
            & (depth_data.position < int(end))
        ]

        region_data ={
                    'primer':primer, 
                    'chromosome':chromosome,
                    'start':start, 
                    'end': end,
                    'depth_df': region_df
        }
        
        data.append(region_data)
    
    return data

def generate_figure(region_data: dict, sample_name: str, output_path: pathlib.Path) -> None: 
    """
    Write a figure that contains all of the depth plots for each target loci. 

    Parameters
    ----------
    region data : dict
        blah
    sample_name : str
        Name of the sample that the output file will have.
    output_path : pathlib.Path
        Path to output figure to.
    """

    def _get_fig_size(num_primers: int) -> int: 
        fig_height = ceil(num_primers/3)
        return fig_height, 3

    num_primers = len(region_data)

    height, width = _get_fig_size(num_primers)

    fig, axs = pyplot.subplots(height, width, figsize=(8.5, 11))

    region_index = 0
    #Iterate through each of the plots
    #Each plot is specified as axs[i][j]
    for i in range(height): 
        for j in range(width): 
            #Only try to generate plots when there's still data
            if region_index < num_primers: 
                
                if region_data[region_index]['depth_df'].empty:
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
                else:
        
                    depth_df = region_data[region_index]['depth_df']
                    
                    print(depth_df.columns)

                    print(depth_df)
                    data = depth_df.loc[:,'count']

                    #Generate bar graph
                    axs[i][j].bar(
                        range(len(depth_df.index)),
                        data,
                        width=1,
                        align='edge',
                    )

                primer = region_data[region_index]['primer']
                axs[i][j].set(
                    title=f'{primer}'
                )
            else: 
                print('NO DATA')
            region_index = region_index + 1
    
    fig_output_path = output_path.joinpath(f'{sample_name}.png')
    fig.subplots_adjust(top=0.94)
    fig.suptitle(sample_name)

    pyplot.tight_layout()

    fig.savefig(fig_output_path, transparent=False, dpi=300)

def parse_args(): 
    parser = argparse.ArgumentParser("Convert tab-delimited values to plots at each target site for read depth.")
    parser.add_argument(
        "input_path",
        action='store', 
        type=pathlib.Path,
        help='Path to tab-separated output of samtools depth for a BAM file.'
    )
    parser.add_argument(
        "primer_bed_path", 
        action='store',
        type=pathlib.Path, 
        help='Path to the primer BED file'
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

    args = parser.parse_args()

    if not(args.output_path): 
        args.output_path = args.input_path
    
    return args

def main(): 
    args = parse_args()

    primer_data = import_primer_data(args.primer_bed_path)

    if args.input_path.is_dir(): 
        for sample_path in args.input_path.glob('*.tsv'):
            sample_name = args.input_path.stem.split('_')[0]

            depth_data = import_depth_data(sample_path)
            region_data = get_region_data(primer_data, depth_data)
            generate_figure(region_data, sample_name, args.output_path)
            
    elif args.input_path.is_file(): 
        sample_name = args.input_path.stem.split('_')[0]

        depth_data = import_depth_data(args.input_path)
        region_data = get_region_data(primer_data, depth_data)
        generate_figure(region_data, sample_name, args.output_path.parent)

    else: 
        print('No tab-separated value files detected.')
        exit()

if __name__ == '__main__': 
    main()