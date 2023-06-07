# core.py

import argparse

def main():
    parser = argparse.ArgumentParser(description='pysift: a Python-based tool for data processing and machine learning')

    parser.add_argument('--input', '-i', type=str, required=True, help='Path to the input file')
    parser.add_argument('--output', '-o', type=str, required=True, help='Path to the output file')
    parser.add_argument('--verbose', '-v', action='store_true', help='If set, increase output verbosity')

    args = parser.parse_args()

    # Call the main logic of your tool using the arguments
    # process_data(args.input, args.output, model=args.model, verbose=args.verbose)

    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    print(f"Verbose: {args.verbose}")

if __name__ == "__main__":
    main()
