import argparse
import csv
import pathlib
import json
import pprint

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
    parser.add_argument(
        '--output_by_sample',
        dest='job_flag',
        action='store_true',
        help='If set, individual json and csvs will be prepared'
    )
    args = parser.parse_args()
    if args.output_path is None: 
        args.output_path = args.json_directory
    return (args.json_directory, args.output_path, args.filter_flag, args.job_flag)

def csv_output(data, output_path): 
    header=(
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
    #Create and setup output CSV file
    csv_output_path = output_path.joinpath('variants.csv')
    csv_file = open(csv_output_path, 'w', newline='')
    csv_writer = csv.writer(csv_file, delimiter=',')
    csv_writer.writerow(header)
    for key in data: 
        csv_writer.writerows(data[key]['variants'])
    csv_file.close()

def json_output(data, output_path): 
    # write to screened JSON file
    prettied_json =  pprint.pformat(data, compact=True, sort_dicts=False, width=160).replace("'", '"')
    with open(output_path.joinpath(f'output_json.json'), 'w') as new_JSON:
        new_JSON.write(prettied_json)

def output_by_loci(json_directory, filter_flag): 
    variant_data = {
    }
    #Generate variants list
    for json_path in json_directory.glob('*.json'):
        json_file = open(json_path, 'r')
        sample_data = json.load(json_file)
        #HARDCODED: CURRENTLY SAMPLE_GENE-PRIMER-DIRECTION ENFORCED IN JSON NAME
        #Distinction between ID and name --> ID does not preserve directionality
        #This is important because we want to pool the SNP calls in reads of forward and reverse direction
        #TO-DO: handle when the snp calls are different between forward and reverse reads
        sample_name = json_path.stem.split('_')[0]
        primer_name = json_path.stem.split('_')[1]
        primer_id = '-'.join(primer_name.split('-')[0:2])
        #Create appropriate data structures
        if primer_id not in variant_data.keys(): 
            variant_data[primer_id] = {'sample_matrix':[], 'variants':[]}
        if sample_name not in variant_data[primer_id]['sample_matrix']: 
            variant_data[primer_id]['sample_matrix'].append(sample_name)
        #Add variant data
        for variant in sample_data['variants']['rows']:
            #[6] - FILTER
            #[2] - ID - populate with gene-primer
            if (
                variant[6] in ('PASS','MANUAL')
                or filter_flag is True
            ):
                variant[2] = sample_name + '_' + primer_name
                variant_data[primer_id]['variants'].append(variant)
        json_file.close()
    return variant_data

def output_by_sample(json_directory, filter_flag): 
    individual_data = {
    }
    #Generate variants list
    for json_path in json_directory.glob('*.json'):
        json_file = open(json_path, 'r')
        sample_data = json.load(json_file)
        #HARDCODED: CURRENTLY SAMPLE_GENE-PRIMER-DIRECTION ENFORCED IN JSON NAME
        #Distinction between ID and name --> ID does not preserve directionality
        #This is important because we want to pool the SNP calls in reads of forward and reverse direction
        #TO-DO: handle when the snp calls are different between forward and reverse reads
        sample_name = json_path.stem.split('_')[0]
        primer_name = json_path.stem.split('_')[1]
        primer_id = '-'.join(primer_name.split('-')[0:2])
        #Create appropriate data structures
        if sample_name not in individual_data.keys():
            individual_data[sample_name] = {'primer_matrix':[], 'variants':[]}
        if primer_name not in individual_data[sample_name]['primer_matrix']: 
            individual_data[sample_name]['primer_matrix'].append(primer_name)
        #Add variant data
        for variant in sample_data['variants']['rows']:
            #[6] - FILTER
            #[2] - ID - populate with gene-primer
            if (
                variant[6] in ('PASS','MANUAL')
                or filter_flag is True
            ):
                variant[2] = sample_name + '_' + primer_name
                individual_data[sample_name]['variants'].append(variant)
        json_file.close()
    return individual_data
        
def main(): 
    json_directory, output_path, filter_flag, job_flag = parse_args()
    if job_flag: 
        variant_data = output_by_sample(json_directory, filter_flag)
    else: 
        variant_data = output_by_loci(json_directory, filter_flag)
    #Generate json output
    json_output(variant_data, output_path)
    #Output csv
    csv_output(variant_data, output_path)

if __name__ == '__main__':
    main()