import argparse
import pathlib
import pandas as pd
import csv
from collections import namedtuple

class CoverageHandler(): 
    #TODO: FILL OUT THIS CLASS WITH THE FUNCTION IN MAIN()
    def __init__(self, coverage_path, primer_path): 
        pass

def parse_args(): 
    parser = argparse.ArgumentParser("Script to analyze the coverage files")
    parser.add_argument(
        'coverage_path',
        action='store', 
        type=pathlib.Path,
        help='Path to coverage files'
    )
    parser.add_argument(
        'primer_path',
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
        args.output_path = args.coverage_path.parent
    return (args.coverage_path, args.primer_path, args.output_path)

def main(): 
    #BIG NOTE:
    #THE PRIMER FILE HAS COORDINATES IN 0-BASED, HALF-OPEN. 
    #THE COVERAGE REPORT IS IN 1-BASED, FULLY-CLOSED.
    #THEREFORE, TO MAINTAIN CONSISTENCY IN MATH, WE WILL CONVERT THE COVERAGE REPORT
    #COORDINATES TO 0-BASED AS WELL. 
    coverage_path, primer_path, output_path = parse_args()
    coverage_df = pd.read_table(coverage_path, names=('chromosome','start','end','methylation_percentage','count_methylated','count_unmethylated'), dtype={'chromosome':'str'})

    #Import primers
    Primer = namedtuple('Primer','chromosome, start, end, primer, sequence, num_cg')
    with open(primer_path, 'r', encoding='utf-8-sig') as primer_file: 
        csv_reader = csv.reader(primer_file, delimiter='\t')
        primer_data = []
        #for primer in map(Primer._make, csv_reader): 
        #    primer_data.append(primer)
        for primer in csv_reader: 
            num_cg = primer[4].count('CG')
            primer_data.append(Primer(
                primer[0],
                primer[1],
                primer[2],
                primer[3],
                primer[4],
                num_cg
            ))
            pass

    #MUTE chained assignment warning
    pd.options.mode.chained_assignment = None
    
    #Summary data
    summary_data = []

    #Perform for each primer set
    for primer in primer_data:
        #Get Dataframe of CG positions within primer region 
        region_df = coverage_df[(coverage_df.chromosome == primer.chromosome) & (coverage_df.start > int(primer.start)) & (coverage_df.end < int(primer.end))]
        #For each CG position, get the position relative to the sequence file in 0-based coordinate, then store in dataframe
        region_df.loc[:,'seq_pos'] = region_df.loc[:,'start'] - 1 - int(primer.start)
        #Sum unmethylated and methylated read calls to get total coverage
        region_df.loc[:,'coverage'] = region_df.loc[:,'count_methylated']+region_df.loc[:,'count_unmethylated']
        #Determine basecall for each CG position (is it C or G)
        basecalls = []
        for cg in region_df.iterrows(): 
            basecalls.append(primer.sequence[cg[1]['seq_pos']])
        region_df.loc[:,'basecall'] = basecalls
        #Check if the current primer set is on the positive or negative strand
        strand = True
        if primer.primer.split('-')[1][0:3] == 'BSN':  
            strand = False
        if strand is True: 
            primer_df = region_df.loc[region_df['basecall']=='C']
        else: 
            primer_df= region_df.loc[region_df['basecall']=='G']
        
        number_cg = len(primer_df)
        methylation_avg = primer_df['methylation_percentage'].mean()
        methylation_std = primer_df['methylation_percentage'].std()
        coverage_avg = primer_df['coverage'].mean()
        coverage_std = primer_df['coverage'].std()
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
        print(f'{primer.primer} mean methylation: {methylation_avg}')
        print(f'{primer.primer} std methylation: {methylation_std}')
        print(f'{primer.primer} mean coverage: {coverage_avg}')
        print(f'{primer.primer} std coverage: {coverage_std}')
        print(f'{primer.primer} number CG under mean coverage: {under_coverage}')

        if primer_df.empty: 
            print('No data.')
        else: 
            print(primer_df)
            primer_data_path = output_path.joinpath(f'{coverage_path.stem}_{primer.primer}.processed.csv')
            primer_df.to_csv(primer_data_path, index=False)

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

    output_csv_path = output_path.joinpath(f'{coverage_path.stem}_summary.csv')
    '''
    with open(output_csv_path, 'w', newline='') as output_csv: 
        csvwriter = csv.writer(output_csv, delimiter=',')
        csvwriter.writerow(
            (
                'primer',
                '#CG_total,
                '#CG_covered',
                'mean_methylation',
                'std_methylation',
                'mean_coverage',
                'std_coverage',
                '#CG<mean',
            )
        )
        csvwriter.writerows(summary_data)
    '''
    summary_df.to_csv(output_csv_path, index=False)

if __name__ == '__main__': 
    main()