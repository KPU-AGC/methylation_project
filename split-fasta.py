import argparse
import pathlib
from Bio import SeqIO

def parse_args(): 
    parser = argparse.ArgumentParser("Split FASTA into several files")
    parser.add_argument(
        "input_path",
        action='store', 
        type=pathlib.Path,
        help='Path to the directory containing input files'
    )
    parser.add_argument(
        '--output',
        dest='output_path',
        action='store',
        default=None,
        type=pathlib.Path,
        help='Output path'
    )

    args = parser.parse_args()
    if args.output_path is None: 
        args.output_path = args.input_path.parent
    return (args.input_path, args.output_path)

def main(): 
    input_path, output_path = parse_args()
    fasta_data = SeqIO.parse(input_path, "fasta")
    for fasta in fasta_data: 
        region = fasta.id.split(':')[0]
        fasta_path = output_path.joinpath(f'{region}.fasta')
        with open(fasta_path, 'w') as fasta_file: 
            SeqIO.write(fasta, fasta_path, 'fasta')

if __name__ == '__main__': 
    main()