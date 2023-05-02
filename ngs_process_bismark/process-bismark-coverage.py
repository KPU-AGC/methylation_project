'''
process-bismark-coverage.py - script to convert coverage files into organized human-readable .csv files

Ensure that the coverage files being provided to this script use 0-based coordinates. 
All calculations performed use 0-based coordinates.

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
import re
from math import ceil
from statistics import mean
#Third-party module
from matplotlib import pyplot
import pandas as pd
import numpy as np

class MethylCallData(): 
    def __init__(self, coverage_path: pathlib.Path, sample_id: str=None): 
        if not(sample_id): 
            self.sample_id = self._determine_id(coverage_path)
        else: 
            self.sample_id = sample_id
        self.methylation_call_data = self._import_coverage_data(coverage_path)
        self.region_data = None
        self.off_target_data = None
        self.summary_data = None

    def _determine_id(self, coverage_path: pathlib.Path) -> str: 
        """Determine the sample ID of this object."""
        return coverage_path.stem.split('_')[0]

    def _import_coverage_data(self, coverage_path: pathlib.Path) -> pd.DataFrame: 
        """Import coverage data into a DataFrame."""
    
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

    def process_methylation_call_data(self, primer_data: dict, return_df: bool=False) -> dict: 
        """
        Generate all of the data structures necessary for downstream analysis

        Parameters
        ----------
        primer_data : list 
            List of Primer namedtuples. 

        Returns
        -------
        primer_dfs : dict
            Key is the primer set, retuns data from target region of primer set.
        """
        def _get_summary_data(df: dict) -> list: 
            """Generate summary data for each primer set."""
            
            summary_data = []

            #Loop through each primer set with data
            for primer in df:
                #Summary statistics 
                num_called_pos = len(df[primer])
                methylation_avg = df[primer]['methylation_percentage'].mean()
                methylation_std = df[primer]['methylation_percentage'].std()
                depth_avg = df[primer]['depth'].mean()
                depth_std = df[primer]['depth'].std()
                under_depth_avg = sum(df[primer]['depth'] < depth_avg)

                #Handling singletons
                if np.isnan(methylation_std): 
                    methylation_std = 0
                if np.isnan(depth_std):
                    depth_std = 0

                summary_data.append(
                    (
                        primer,
                        num_called_pos, 
                        methylation_avg, 
                        methylation_std, 
                        depth_avg,
                        depth_std,
                        under_depth_avg,
                    )
                )

            return summary_data

        def _get_off_target_df(primers: dict, df: pd.DataFrame) -> pd.DataFrame: 
            '''Figure out what other sequences remain.'''
            off_target_df = df
            #print(off_target_df)
            for primer in primers.values(): 
                off_target_df.drop(
                    off_target_df[
                        (off_target_df.chromosome == primer.chromosome)
                        & (off_target_df.start >= int(primer.start))
                        & (off_target_df.start <= int(primer.end))
                    ].index,
                    inplace=True,
                )

            return off_target_df

        primer_dfs = dict()

        methylation_df = self.methylation_call_data

        #Perform for each primer set
        for primer in primer_data.values():

            #Get Dataframe of methylation called positions within primer region 
            region_df = methylation_df[
                (methylation_df.chromosome == primer.chromosome) 
                & (methylation_df.start >= int(primer.start)) 
                & (methylation_df.end <= int(primer.end))
            ]

            if region_df.empty:
                print(f'No data for {primer.primer}')
            else:
                #For each CG position, get the position relative to the sequence file 
                # in 0-based coordinate, then store in dataframe.
                #Assumes input coverage file was 0-based.
                region_df['seq_pos'] = region_df['start'] - int(primer.start)

                #Sum unmethylated and methylated read calls to get total coverage
                region_df['depth'] = (
                    region_df['count_methylated'] 
                    + region_df['count_unmethylated']
                )

                #Determine basecall for each called position (C or G)
                basecalls = []
                seq_pos = region_df.loc[:, 'seq_pos'].tolist()
                for pos in seq_pos: 
                    basecalls.append(primer.sequence[pos])
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

                primer_dfs[primer.primer] = primer_df

        self.region_data = primer_dfs
        self.off_target_data = _get_off_target_df(primer_data, self.methylation_call_data)
        self.summary_data = _get_summary_data(primer_dfs)

        if return_df: 
            return primer_dfs

    def generate_figure(
        self, 
        primer_data: list, 
        output_path: pathlib.Path, 
        cg_flag: bool=False,
        show_fig_flag: bool=False,
    ): 
        """
        Generate a figure containing methylation plots for each target site. 

        Parameters
        ----------
        primer_data : dict
            key (primer set) yields information about the primer set
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

        def _prep_region_data(region_data: dict, cg_flag: bool, primer_data: dict=None) -> list: 
            """
            Convert the DataFrame data into a form that is easier to plot in matplotlib.

            Parameters
            ----------
            region_data : dict
            cg_flag : bool
            primer_data : dict

            Returns
            -------
            target_sites : list
                List of dictionaries where each dictionary is one target site, and has the plotted
                information.

            """

            target_sites = []
            for key in region_data: 
                target_data = {
                    'primer':key,
                    'positions':region_data[key].loc[:,'start'].to_list(),
                    'methylation':region_data[key].loc[:,'methylation_percentage'].to_list(),
                    'depth':region_data[key].loc[:,'depth'].to_list(),
                }
                if cg_flag:
                    #Include only CG positions if 
                    cg_pos = primer_data[key].cg_pos

                    only_cg_pos = []
                    only_methyl_pos = []
                    only_depth_pos = []

                    #Check each position if it's a CG.
                    #Append data if it is. 
                    for index, pos in enumerate(target_data['positions']): 
                        #Conditional allows checking of BPX and BNX primers (C or G)
                        if (pos in cg_pos) or (pos+1 in cg_pos): 
                            only_cg_pos.append(pos)
                            only_methyl_pos.append(target_data['methylation'][index])
                            only_depth_pos.append(target_data['depth'][index])
                        else:
                            print(f'{pos} not found.')

                    target_data['positions'] = only_cg_pos
                    target_data['methylation'] = only_methyl_pos
                    target_data['depth'] = only_depth_pos
                try: 
                    #Calculating relative coverage
                    max_depth = max(target_data['depth'])
                    target_data['relative_depth'] = [(depth/max_depth)*100 for depth in target_data['depth']]
                except ValueError: 
                    print('No data')
                target_sites.append(target_data)

            return target_sites

        region_data = _prep_region_data(self.region_data, cg_flag, primer_data)

        num_primers= len(primer_data)
        height, width = _get_fig_size(num_primers)
        fig, axs = pyplot.subplots(height, width, figsize=(8.5,11))

        region_index = 0

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
                        
                        max_depth = max(region_data[region_index]['depth'])

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
                            region_data[region_index]['relative_depth'],
                            'red'
                            )
                        axs[i][j].text(
                            0,
                            3,
                            f'max(depth)={str(max_depth)}',
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
        fig_output_path = output_path.joinpath(f'{self.sample_id}.png')
        fig.subplots_adjust(top=0.94)
        fig.suptitle(self.sample_id)
        #pyplot.tight_layout()

        fig.savefig(fig_output_path, transparent=False, dpi=300)

        if show_fig_flag: 
            pyplot.show()

    def export_csv(self, primer_path: pathlib.Path, summary_path: pathlib.Path) -> None: 
        """"""
        
        #Output primer DataFrames
        for primer in self.region_data: 
            primer_data_path = primer_path.joinpath(f'{self.sample_id}_{primer}_processed.csv')
            self.region_data[primer].to_csv(primer_data_path, index=False)
        #Output summary data
        summary_names=(
            'primer',
            '#covered_positions',
            'mean_methylation',
            'std_methylation',
            'mean_depth',
            'std_depth',
            '#CG<mean',
        )
        
        output_csv_path = summary_path.joinpath(f'{self.sample_id}_summary.csv')
        with open(output_csv_path, 'w', newline='',) as csv_file: 
            csv_writer = csv.writer(csv_file)
            csv_writer.writerow(summary_names)
            csv_writer.writerows(self.summary_data)

def import_primer_data(path: pathlib.Path) -> dict(): 
    """Return primer data from bed file specifiying target regions"""

    primer_data = dict()

    Primer = namedtuple(
        'Primer',
        'chromosome, start, end, primer, sequence, num_cg, cg_pos'
    )

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
                #1-based index - should only have the C positions
                cg_locations.append(cg.start()+int(primer[1])+1)
            
                primer_data[primer[3]] = Primer(
                    primer[0],
                    primer[1],
                    primer[2],
                    primer[3],
                    primer[4],
                    num_cg, 
                    cg_locations,
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
    parser.add_argument(
        '-s', 
        '--snpsplit',
        dest='snp_flag',
        action='store_true',
        help='Set this parameter if you are processing several bams per file (SNPsplit)',
    )
    parser.add_argument(
        '-c',
        '--cg_only',
        dest='cg_flag',
        action='store_true',
        help='Set this parameter if you want the charts to only show CG methylation.'
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

    if args.coverage_path.is_dir(): 
        coverage_paths = get_coverage_paths(args.coverage_path)

        summary_dir = args.output_path.joinpath('summaries')
        primer_dir = args.output_path.joinpath('individual_primer')
        figure_dir = args.output_path.joinpath('summary_fig')
        summary_dir.mkdir(exist_ok=True)
        primer_dir.mkdir(exist_ok=True)
        figure_dir.mkdir(exist_ok=True)

        for coverage_path in coverage_paths: 
            if args.snp_flag: 
                sample_name = coverage_path.stem.split('_')[0]
                #Check for genome assignment
                if 'genome1' in coverage_path.stem: 
                    pool_name = 'genome1'
                elif 'genome2' in coverage_path.stem: 
                    pool_name = 'genome2'
                elif 'unassigned' in coverage_path.stem:
                    pool_name = 'unassigned'
                else: 
                    pool_name ='all'
                #pool_name = coverage_path.stem.split('.')[-2]
                sample_id = f'{sample_name}_{pool_name}'
                mcd = MethylCallData(coverage_path, sample_id)
            else: 
                mcd = MethylCallData(coverage_path)
            
            mcd.process_methylation_call_data(primer_data)
            mcd.generate_figure(primer_data, figure_dir, args.cg_flag)
            mcd.export_csv(primer_dir, summary_dir)
      
    elif args.coverage_path.is_file(): 
        primer_dir = args.output_path.joinpath('individual_primer')
        primer_dir.mkdir(exist_ok=True)

        mcd = MethylCallData(args.coverage_path)
        mcd.process_methylation_call_data(primer_data)
        mcd.generate_figure(primer_data, args.output_path, args.cg_flag)
        mcd.export_csv(primer_dir, args.output_path)

    else:
        print('Path is neither a directory or file.')
        exit()

if __name__ == '__main__': 
    main()