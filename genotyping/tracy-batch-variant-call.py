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
import csv
from Bio import SeqIO

class Args(NamedTuple):
    """ Command-line arguments """
    reference_file_path: pathlib.Path
    target_path: pathlib.Path
    pratio: float
    output: pathlib.Path

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
        dest='output',
        type=pathlib.Path,
        default=None,
        metavar='path',
        help='the path of the output')

    args = parser.parse_args()

    args.target_path = pathlib.Path.resolve(args.target_path)

    if not args.output: 
        args.output = args.target_path.joinpath('output')
        print(args.output)
        args.output.mkdir(exist_ok=True)
    
    if args.pratio < 0 or args.pratio > 1:
        parser.error('Peak ratio must be between 0 and 1.')

    if not args.target_path.is_dir():
        parser.error('The input must be a directory!')

    return Args(args.reference_file_path, args.target_path, args.pratio, args.output)


def run_tracy(primer_id_arg: str, sample_id_arg: str, reference_fasta_path_arg: pathlib.Path, ab1_file_path_arg: pathlib.Path, peak_ratio_arg: float, trim_left: int, trim_right: int):
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
        '0',
        '-q',
        f'{trim_left}',
        '-u',
        f'{trim_right}',
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
        input_path = target_path_arg.joinpath(f'{input_prefix_arg}{suffix}')
        output_path = output_path_arg.joinpath(f'{input_prefix_arg}{suffix}')
        try:
            shutil.move(input_path, output_path)
        except FileNotFoundError:
            pass
        finally:
            print(input_path)
            print(output_path)


def check_name(path):
    """
    Checks if the name is ok
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


def abi_trim(ab1_file_path):
    untrimmed = SeqIO.read(ab1_file_path, 'abi')
    trimmed = SeqIO.read(ab1_file_path, 'abi-trim')
    left_trim = untrimmed.seq.find(trimmed.seq[0:5])-1
    right_trim = len(untrimmed.seq) - len(trimmed) - left_trim
    return left_trim, right_trim


def output_summary(output_path, best_runs, failed_runs): 
    csv_output_path = output_path.joinpath("ab1_summary.csv")
    csv_file = open(csv_output_path, 'w')
    csv_writer = open(csv_file, delimiter=',')
    csv_writer.writerow('Best runs')
    csv_writer.writerow(('Sample', 'Trace score', 'Median PUP'))
    csv_writer.writerows(best_runs)
    csv_writer.writerow('Failed runs')
    csv_writer.writerow(('Sample', 'Trace score', 'Median PUP'))
    csv_writer.writerows(failed_runs)
    csv_file.close()


# --------------------------------------------------
def main() -> None:
    """ Perform the function """

    args = get_args()
    target_path = pathlib.Path.resolve(args.target_path)

    #output_path = target_path.parent.joinpath('tracy-files')
    #output_path.mkdir(parents=True, exist_ok=True)

    list_of_samples = {}
    failed_samples = []
    best_run_samples = []

    #Organize the ab1 files
    #List of unique groups containing sequences for a particular sample and primer
    for sample in target_path.glob('*.ab1'):
        if '_'.join(list(check_name(sample).values())) not in list_of_samples:
            list_of_samples['_'.join(list(check_name(sample).values()))] = []

        scores = get_ab1_qc(sample)

        if scores['trace_score'] >= 30 and scores['median_PUP_score'] >= 10:
            list_of_samples['_'.join(list(check_name(sample).values()))].append(sample)
        else: 
            failed_samples.append((sample.name, scores['trace_score'], scores['median_PUP_score']))

    #Identify the sequence with the best trace score
    for sample in list_of_samples:
        best_run = None
        best_trace = 0
        for run in list_of_samples[sample]:
            scores = get_ab1_qc(run)
            if scores['trace_score'] > best_trace:
                best_run = run
                best_trace = scores['trace_score']
                best_pup = scores['median_PUP_score']

        if best_run:
            run_name = check_name(best_run)
            left_trim, right_trim = abi_trim(best_run)
            best_run_samples.append((run.name, best_trace, best_pup))

            run_tracy(
                primer_id_arg=f"{run_name['primer_id']}-{run_name['direction']}",
                sample_id_arg=run_name['sample_id'],
                reference_fasta_path_arg=args.reference_file_path,
                ab1_file_path_arg=best_run,
                peak_ratio_arg=args.pratio,
                trim_left=left_trim,
                trim_right=right_trim)

            if (args.output.stem == 'output'):
                move_tracy_files(pathlib.Path.cwd(), f"{run_name['sample_id']}_{run_name['primer_id']}-{run_name['direction']}", pathlib.Path(args.output))
                
    output_summary(args.output, best_run_samples, failed_samples)
# --------------------------------------------------
if __name__ == '__main__':
    main()
