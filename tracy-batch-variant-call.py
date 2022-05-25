#!/usr/bin/env python3
"""
Author : Erick Samera
Date   : 2022-05-22
Purpose: To batch together 'tracy decompose' commands.
"""

import argparse
from typing import NamedTuple, TextIO
import time

import subprocess
import pathlib
import shutil
from Bio import SeqIO

class Args(NamedTuple):
    """ Command-line arguments """
    reference_file_path: pathlib.Path
    target_path: pathlib.Path
    pratio: float
    trim: int
    output_path: pathlib.Path

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description="Perform 'tracy assemble' on a group of files.",
        epilog='',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        'reference_file_path',
        type=pathlib.Path,
        metavar='path',
        help='path to reference fasta/genome')

    parser.add_argument(
        'target_path',
        type=pathlib.Path,
        metavar='path',
        help='path to directory containing .ab1 files')
    
    parser.add_argument(
        '-p',
        '--pratio',
        type=float,
        metavar='<n>',
        default=0.8,
        help='tracy peak ratio, [0-1]')

    parser.add_argument(
        '-o',
        '--output',
        type=pathlib.Path,
        metavar='path',
        help='the path of the output')

    parser.add_argument(
        '-t',
        '--trim',
        type=int,
        metavar='<n>',
        default=2,
        help='tracy trim stringency {1:9}, 1 = most')

    args = parser.parse_args()

    args.target_path = pathlib.Path.resolve(args.target_path)
    
    if args.pratio < 0 or args.pratio > 1:
        parser.error('Peak ratio must be between 0 and 1.')

    if args.trim not in range(1,10):
        parser.error('Trim stringency value must be int between 1 and 9, inclusive.')

    if not args.target_path.is_dir():
        parser.error('The input must be a directory!')

    return Args(args.reference_file_path, args.target_path, args.pratio, args.trim, args.output_path)


def run_tracy(primer_id_arg: str, sample_id_arg: str, reference_fasta_path_arg: pathlib.Path, ab1_file_path_arg: pathlib.Path, peak_ratio_arg: float, trim_arg: int):
    """
    Run tracy using the reference fasta files
    """

    args=[
        'tracy',
        'decompose',
        '-r',
        f'{reference_fasta_path_arg}',
        '-v',
        '-t',
        f'{trim_arg}',
        '-p',
        f'{peak_ratio_arg}',
        '-o',
        f'{sample_id_arg}_{primer_id_arg}',
        ab1_file_path_arg,
    ]
    result = subprocess.run(args, capture_output=False)


def move_tracy_files(target_path_arg: pathlib.Path, input_prefix_arg: str, output_path_arg: pathlib.Path) -> None:
    """
    Move files produced by 'tracy assemble' into another folder for easier file management.

    Parameters:
        target_path_arg (pathlib.Path): path containing .ab1 files
        input_prefix_arg (str): sample name prefix for files
        output_path_arg (pathlib.Path): directory containing

    Returns:
        None
    """

    # suffix list of files produced by tracy assemble
    tracy_files_suffix_list = ['.abif', '.align1', '.align2', '.align3', '.bcf', '.bcf.csi', '.decomp', '.json']

    # move each file into the output folder
    for suffix in tracy_files_suffix_list:
        try:
            shutil.move(target_path_arg.joinpath(f'{input_prefix_arg}{suffix}'), output_path_arg.joinpath(f'{input_prefix_arg}{suffix}'))
        except FileNotFoundError:
            pass


def check_name(path):
    """
    asdfasfd
    """

    if not isinstance(path, pathlib.Path): 
        raise TypeError("Input is not a path")

    name_split = path.stem.split('_')
    #Determine which filename convention it is
    index = 0
    for item in name_split:
        if len(item) == 2:
            index = name_split.index(item)
            break
        print(item)
    #Check if there is primer direction specified
    primer_info = name_split[1].split('-')
    if len(primer_info) == 3:
        primer_direction = primer_info[2]
    else: 
        primer_direction = 'F'
    #Store components of filename
    sequence_info = {
        0:{'sample_id':name_split[2], 'primer_id':'-'.join(primer_info[0:2]),'direction':primer_direction},
        1:{''},
        2:{'sample_id':name_split[0], 'primer_id':'-'.join(primer_info[0:2]),'direction':primer_direction},
    }
    return sequence_info[index]


def get_ab1_qc(path_arg: pathlib.Path) -> dict:
    """
    Get the QC data for each ab1 file
    """

    query_ab1_file = SeqIO.read(path_arg, 'abi')

    try:
        return {'trace_score': int(query_ab1_file.annotations['abif_raw']['TrSc1']), 'median_PUP_score': int(query_ab1_file.annotations['abif_raw']['PuSc1'])}
    except:
        return {'trace_score': 0, 'median_PUP_score': 0}


def print_runtime(action) -> None:
    """Prints the runtime."""
    print(f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}] {action}')

# --------------------------------------------------
def main() -> None:
    """ Perform the function """

    args = get_args()
    target_path = pathlib.Path.resolve(args.target_path)

    #output_path = target_path.parent.joinpath('tracy-files')
    #output_path.mkdir(parents=True, exist_ok=True)

    list_of_samples = {}

    for sample in target_path.glob('*.ab1'):
        if '_'.join(list(check_name(sample).values())) not in list_of_samples:
            list_of_samples['_'.join(list(check_name(sample).values()))] = []
        if get_ab1_qc(sample)['trace_score'] >= 30 and get_ab1_qc(sample)['median_PUP_score'] >= 10:
            list_of_samples['_'.join(list(check_name(sample).values()))].append(sample)

    for sample in list_of_samples:
        best_run = None
        best_trace = 0
        for run in list_of_samples[sample]:
            if get_ab1_qc(run)['trace_score'] > best_trace:
                best_run = run
                best_trace = get_ab1_qc(run)['trace_score']

        if best_run:
            run_data = check_name(best_run)

            run_tracy(
                primer_id_arg=f"{run_data['primer_id']}-{run_data['direction']}",
                sample_id_arg=run_data['sample_id'],
                reference_fasta_path_arg=args.reference_file_path,
                ab1_file_path_arg=best_run,
                peak_ratio_arg=args.pratio,
                trim_arg=args.trim)

            move_tracy_files(best_run.parent, f"{run_data['sample_id']}_{run_data['primer_id']}-{run_data['direction']}", args.output_path)

# --------------------------------------------------
if __name__ == '__main__':
    main()
