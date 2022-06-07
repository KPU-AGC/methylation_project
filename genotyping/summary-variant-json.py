import argparse
import csv
import pathlib
import json
import pprint
import pandas as pd

def parse_args(): 
    parser = argparse.ArgumentParser('Output the variant data for each site in a human-readable format')
    parser.add_argument(
        'json_path',
        action='store',
        type=pathlib.Path,
        help='Path to the json output of combine-variant-calls.py'
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
        '--template',
        dest='template_path',
        action='store',
        default=None,
        type=pathlib.Path,
        help='Path to template file',
    )
    args = parser.parse_args()
    if args.output_path is None: 
        args.output_path = args.json_path.parent
    return args.json_path, args.output_path, args.template_path

def json_output(data, output_path): 
    # write to screened JSON file
    prettied_json =  pprint.pformat(data, compact=True, sort_dicts=False, width=160).replace("'", '"')
    with open(output_path.joinpath(f'human_table.json'), 'w') as new_JSON:
        new_JSON.write(prettied_json)

def get_template_index(template_path):
    with open(template_path, 'r') as template_file: 
        return template_file.read().split('\n')

def main(): 
    json_path, output_path, template_path = parse_args()

    json_file = open(json_path, 'r')
    sample_data = json.load(json_file)
    json_file.close()

    all_data = dict()
    #Process all variants
    for gene in sample_data:
        variant_data = dict()
        for variant in sample_data[gene]['variants']:
            #Generate chromosome position
            variant_position = 'chr'+variant[0]+':'+str(variant[1])
            sample_name = variant[2].split('_')[0]
            #Add data to variant_data
            if variant_position not in variant_data: 
                #set reference; assumes haploid reference
                variant_data[variant_position] = {
                    'REF':variant[3]+'/'+variant[3],
                }
            #Check variant type
            if variant[8] == 'hom. ALT':
                variant_data[variant_position][sample_name] = variant[4]+'/'+variant[4]
            elif variant[8] == 'het.':
                variant_data[variant_position][sample_name] = variant[3]+'/'+variant[4]
            else:
                raise ValueError(f"Variant type not hom. ALT or het.\n{variant}")
        #Now that all variants are processed, add samples that match reference to the data table
        for variant_position in variant_data:
            for sample_name in sample_data[gene]['sample_matrix']:
                if sample_name not in variant_data[variant_position].keys(): 
                    variant_data[variant_position][sample_name] = variant_data[variant_position]['REF']
        #Process into dataframe, reindex
        if template_path is None: 
            index_df = sorted(sample_data[gene]['sample_matrix'])
        else: 
            index_df = get_template_index(template_path)
        variant_df = pd.DataFrame(variant_data).reindex(index_df)
        all_data[gene] = variant_df

    #Output file
    for gene in all_data:
        with open(output_path.joinpath(f'{gene}_variants.csv'), 'w', newline='') as output_csv:
            all_data[gene].to_csv(output_csv, na_rep='NaN')

if __name__ == '__main__': 
    main()