#!/usr/bin/env python3
__description__ =\
"""
Purpose: To simulate targeted-amplicon bisulfite-converted reads with given experimental weights.
"""
__author__ = "Erick Samera"
__version__ = "0.0.1"
__comment__ = 'stable enough??'

# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    RawTextHelpFormatter)
import sys
import pathlib
import pickle
import subprocess
import tempfile
import pandas as pd
import random
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description=f"{__description__}",
        epilog=f"v{__version__} : {__author__} | {__comment__}",
        formatter_class=RawTextHelpFormatter)

    req_args = parser.add_argument_group('Required Args (*)')
    req_args.add_argument(
        '-w',
        '--weights-path',
        metavar='<PKL>',
        type=str,
        default=None,
        help="(*) Path to .pkl containing weights per amplicon.")
    req_args.add_argument(
        '-r',
        '--region',
        metavar='<STR>',
        dest='region',
        default=None,
        help="(*) Region name for generating amplicon. (ex. DNMT3A-BP1X)\nUse 'list' to show options.")
    parser.add_argument(
        '-m',
        '--metadata',
        metavar='<PKL>',
        type=pathlib.Path,
        required=False,
        help="Path to .pkl containing metadata. [None]")
    parser.add_argument(
        '-f',
        '--apply-filter',
        metavar='<STR>',
        dest='filter_condition',
        required=False,
        help="Condition to filter for applying weights, case-sensitive (ex. survived_birth=TRUE or large_calf=FALSE) [None]\nUse 'list' to show options.")

    sim_args = parser.add_argument_group('Simulation Args', 'WIP: math not mathing, defaults result in 600X')
    sim_args.add_argument(
        '-s',
        '--subsample',
        metavar='<INT>',
        type=int,
        default=15,
        help="Subsampling depth in the list of allele fractions for weights and amplicon creation. [15]")
    sim_args.add_argument(
        '-S',
        '--seq-depth',
        metavar='<INT>',
        type=int,
        default=100,
        help="Sequencing depth. [100]")
    sim_args.add_argument(
        '-W',
        '--manual-weights',
        metavar='<str>',
        type=str,
        help="Apply manual weights. [None]\nex:\n - \"100\" = all reads fully methylated.\n - \"100:0.5,0:0.5\" = half reads fully methylated, other half not methylated.")

    file_args = parser.add_argument_group('Output Args')
    file_args.add_argument(
        '-p',
        '--out-prefix',
        metavar='<STR>',
        type=str,
        default='SIMULATED',
        help="Sample prefix of output {OUT-PREFIX}.{REGION}.R{1/2}.fastq.gz [SIMULATED]")
    file_args.add_argument(
        '-o',
        '--out-dir',
        metavar='<PATH>',
        type=pathlib.Path,
        default=None,
        help="Output directory. [in-place]")

    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------
    if not args.out_dir: args.out_dir = pathlib.Path().cwd()

    if args.region:
        if args.region.lower() == 'list':
            choices=["DMAP1-BP1X","DNMT3A-BP1X","DNMT3B-BP2X","GNAS-BN1X","H19-BP1X","IGF2R-BN2X","KCNQ1-BP1X","LIF-BN1X","LIFR-BP1X","MEST-BN1X","NNAT-BP1X","PEG10-BP1X","PEG3-BN1X","PLAGL1-BN1X","RTL1-BN1X","SLC2A8-BP1X","SNRPN-BP1X","SUV39H1-BN1X","TXNIP-BP1X","XIST-BP1X"]
            print('Regions list:')
            for i in choices: print('-',i)
            quit()
    if args.filter_condition:
        if args.filter_condition.lower() == 'list':
            choices=[
                "survived_birth:\tTRUE/FALSE",
                "has_some_abnormality:\tTRUE/FALSE",
                "control:\tTRUE/FALSE",
                "large_calf:\tTRUE/FALSE",
                "sex:\tF/M/?"]
            print('Regions list:')
            for i in choices: print('-',i)
            quit()

    if (not args.weights_path) or (not args.region):
        parser.print_help()
        print()
        print('ERROR: Both -w/--weights-path and -r/--region are required!')
        quit()
    return args
# --------------------------------------------------
def _simulate_reads(args, sequence_weights: list, output_dir: pathlib.Path) -> None:
    """
    Simulates reads based on the provided sequence weights.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        fasta_output = pathlib.Path(temp_dir).joinpath('templates')
        fasta_output.mkdir(exist_ok=True, parents=True)

        fastq_output = pathlib.Path(temp_dir).joinpath('reads')
        fastq_output.mkdir(exist_ok=True, parents=True)

        generated_sequences = random.choices([i[0] for i in sequence_weights], [i[1] for i in sequence_weights], k=100)
        for i, seq in enumerate(generated_sequences):
            current_iter = fasta_output.joinpath(f'iter.{i}.fasta')
            with open(current_iter, mode='w') as output_file:
                output_file.write(f'>iter.{i}\n{seq}\n')
            _do_art_simulation(args, current_iter, len(seq), fastq_output)
        _combine_reads(args, fastq_output, output_dir)
    return None
def _do_art_simulation(args, template_file: pathlib.Path, read_len: int, out_dir: pathlib.Path) -> None:
    """
    Simulates sequencing reads using ART.
    """
    iteration: str = template_file.stem.split('.')[1]
    seq_depth: int = int(args.seq_depth / args.subsample)

    #read_len = '-l 250' if read_len > 250 else read_len
    command = f'art_illumina -ss MSv3 -amp -p -na -i {template_file} -m {read_len} -l 250 -f {seq_depth} -o {out_dir.joinpath(f'iter.{iteration}.R')}'
    subprocess.run(command, shell=True, capture_output=True)
    return None
def _combine_reads(args, input_dir: pathlib.Path, out_dir: pathlib.Path) -> None:
    """
    """
    for read_n in ("R1", "R2"):
        files_to_concat = [str(file) for file in input_dir.glob(f'*{read_n}.fq')]
        command = f'cat ' + ' '.join(files_to_concat) + ' | gzip > ' + str(out_dir.joinpath(f'{args.out_prefix}.{args.region}.{read_n}.fastq.gz'))
        subprocess.run(command, shell=True, capture_output=True)
def _aggregate_read_stats(input_normalized_df) -> list:
    """
    """
    column_stats = {}
    for column in input_normalized_df.columns:
        if column != 'sample.id':
            non_zero_values = input_normalized_df[column][input_normalized_df[column] != 0]
            avg_value = non_zero_values.mean()
            std_dev = non_zero_values.std()
            column_stats[column] = {
                'average': avg_value,
                'std_dev': std_dev,
            }

    stats_list = [{'column': col, **stats} for col, stats in column_stats.items()]
    return sorted(stats_list, key=lambda x:x['average'], reverse=False)
def _generate_weights(args, amplicon_weights, metadata=None) -> pd.DataFrame:
    """
    """
    
    if (args.metadata) and (args.filter_condition):
        filter_attr, filter_val = args.filter_condition.split('=')
        weights_filtered = amplicon_weights[args.region][amplicon_weights[args.region]['sample.id'].isin(metadata[filter_attr][filter_val])]
    else:
        weights_filtered = amplicon_weights[args.region]

    normalized_df = weights_filtered.drop(columns=['sample.id']).div(weights_filtered.drop(columns=['sample.id']).sum(axis=1), axis=0)
    normalized_df['sample.id'] = weights_filtered['sample.id']
    normalized_df = normalized_df[['sample.id'] + [col for col in normalized_df.columns if col != 'sample.id']]
    normalized_df.fillna(-1)
    return normalized_df
def _process_amplicon(args, CpG_methylation_p: int, amplicon_weights) -> None:
    """
    """
    template_sequence = amplicon_weights['templates'][args.region]
    CpG_positions = [i for i, nuc in enumerate(template_sequence) if nuc.upper() == 'C']
    num_unmethylated = int(len(CpG_positions) * (1 - CpG_methylation_p / 100))
    unmethylated = sorted(random.sample(CpG_positions, k=num_unmethylated))
    
    output_sequence = list(template_sequence)
    for i in unmethylated:
        output_sequence[i] = 'T'
    output_sequence = ''.join(output_sequence)
    return output_sequence
def _parse_custom_weights(args, weights_str: str, amplicon_weights=None) -> None:
    """
    """

    def _print_invalid(reason_str: str) -> None:
        """
        """
        print('ERROR! Invalid weight string!')
        print(reason_str)
        quit()

    try:
        weights_list = [(int(i.split(':')[0]), float(i.split(':')[1])) if ':' in i else int(i) for i in weights_str.split(',')]
        if ':' in weights_str:
            if not all([0 <= i[0] <= 100 for i in weights_list]): _print_invalid('Methylation percentage must be between 0 and 100.')
            if not sum([i[1]for i in weights_list]) == 1: _print_invalid('Proportion of reads should add up to 100%.')
            sequences = [_process_amplicon(args, i[0], amplicon_weights) for i in weights_list]
            weights = [i[1] for i in weights_list]
        else:
            sequences = [_process_amplicon(args, i, amplicon_weights) for i in weights_list]
            weights = [100]
    except:
        quit()
    
    
    return sequences, weights
def errprint(*args, **kwargs): print(*args, file=sys.stderr, **kwargs)
# --------------------------------------------------
def main() -> None:
    """ Insert docstring here """

    args = get_args()

    try:
        with open(pathlib.Path(args.weights_path), mode='rb') as input_file: amplicon_weights: dict = pickle.load(input_file)
        amplicon_weights.get('KCNQ1-BP1X')
    except KeyError: print('oopsies')

    if args.metadata:
        try:
            with open(args.metadata, mode='rb') as input_file: metadata: dict = pickle.load(input_file)
        except KeyError: print('oopsies')

    if not args.manual_weights:
        normalized_weights: pd.DataFrame = _generate_weights(args, amplicon_weights, metadata) if args.metadata else _generate_weights(args, amplicon_weights)
        read_statistics: list = _aggregate_read_stats(normalized_weights)[:args.subsample]

        weights: list = [random.uniform(i['average']-i['std_dev'], i['average']+i['std_dev']) for i in read_statistics]
        sequences: list = [i['column'] for i in read_statistics]
        sequence_weights: list = [(sequence, weight) for sequence, weight in zip(sequences, weights) if weight>0]
    elif args.manual_weights:
        sequences, weights = _parse_custom_weights(args, args.manual_weights, amplicon_weights)
        sequence_weights: list = [(sequence, weight) for sequence, weight in zip(sequences, weights)]

    errprint(f'Simulating reads from {args.region} to {args.out_dir}/{args.out_prefix}.{args.region}.'+"{R1/R2}"+'.fastq.gz ...')
    _simulate_reads(args, sequence_weights, args.out_dir)
    errprint('Done!')

    return None
# --------------------------------------------------
if __name__ == '__main__': 
    main()
