import argparse
import pathlib
import pandas as pd
import csv
from collections import namedtuple

class CoverageHandler(): 
    #TODO: FILL OUT THIS CLASS WITH THE FUNCTION IN MAIN()
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
    coverage_df = pd.read_table(coverage_path, names=('chromosome','start','end','methylation_percentage','count_methylated','count_unmethylated'))
    #Import primers
    Primer = namedtuple('Primer','chromosome, start, end, primer, sequence')
    with open(primer_path, 'r', encoding='utf-8-sig') as primer_file: 
        csv_reader = csv.reader(primer_file, delimiter='\t')
        primer_data = []
        for primer in map(Primer._make, csv_reader): 
            primer_data.append(primer)
    #MUTE chained assignment warning
    pd.options.mode.chained_assignment = None
    #Perform for each primer set
    for primer in primer_data:
        #Get Dataframe of CG positions within primer region 
        region_df = coverage_df[(coverage_df.chromosome == primer.chromosome) & (coverage_df.start > int(primer.start)) & (coverage_df.end < int(primer.end))]
        #Get position on sequence file, in 0-based coordinates, for each CG position
        region_df.loc[:,'seq_pos'] = region_df.loc[:,'start'] - 1 - int(primer.start)
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
        
        #for cg in primer_df.iterrows(): 
        #    cg_data = cg[1]
            #Coordinate of the CG is 1-based.
            #Convert to 0, then subtract to find the 0-based coordinate, starting at 0. 
        #    cg_pos = cg_data['start'] - 1 - int(primer.start)
            
            
            
        #    start = cg_data['start']
        #    print(f'{start} : {primer.start}')
        #    print(f'{primer.sequence[cg_pos]} : {primer.sequence[0]} : {primer.primer}')
        #print(primer.sequence)
        #print(primer_df)
        methylation_avg = primer_df['methylation_percentage'].mean()
        methylation_std = primer_df['methylation_percentage'].std()
        methylation_count_avg = primer_df['count_methylated'].mean()
        methylation_count_std = primer_df['count_methylated'].std()
        unmethylated_count_avg = primer_df['count_unmethylated'].mean()
        unmethylated_count_std = primer_df['count_unmethylated'].std()
        coverage_avg = (primer_df['count_methylated'] + primer_df['count_unmethylated']).mean()
        coverage_std = (primer_df['count_methylated'] + primer_df['count_unmethylated']).std()
        print(f'{primer.primer} mean methylation: {methylation_avg}')
        print(f'{primer.primer} std methylation: {methylation_std}')
        print(f'{primer.primer} mean coverage: {coverage_avg}')
        print(f'{primer.primer} std coverage: {coverage_std}')

        if primer_df.empty: 
            print('No data.')
        else: 
            print(primer_df)

if __name__ == '__main__': 
    main()