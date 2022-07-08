import argparse
import pathlib

def parse_args(): 
    parser = argparse.ArgumentParser("Plcaeholder")
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
        args.output_path = args.input_path
    return (args.input_path, args.output_path)

def main(): 
    input_path, output_path = parse_args()

if __name__ == '__main__': 
    main()