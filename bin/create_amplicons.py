from ncov.primers import *
import os
import argparse

def init_args():
    """
    Initialize the command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', required=True,
                        help='BED file containing primers to process')
    parser.add_argument('-p', '--path', default='./qc_primers',
                        help='')
    parser.add_argument('-o', '--output', default='amplicon.bed',
                        help='BED file to write amplicons')
    parser.add_argument('-t', '--type', default='unique_amplicons',
                        help='Type of amplicons to generate (unique_amplicons, full, no_primers)')
    return parser.parse_args()


def main():
    """
    Main program
    """
    args = init_args()
    primers = import_bed_file(bed=args.file)
    primer_dict = dict()
    for primer in primers:
        add_primer(primers=primer_dict, primer=primer)
    amplicons = create_amplicon(primers=primer_dict, type=args.type,
                                path=args.path)
    write_amplicon_to_bed(amplicons=amplicons, outfile=args.output)
    print(f"BED file written to {args.output}")


if __name__ == '__main__':
    main()

