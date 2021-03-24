#!/usr/bin/env python

from ncov.primers import *
import os
import argparse
import shutil
import tempfile

def init_args():
    """
    Initialize the command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', required=True,
                        help='BED file containing primers to process')
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
    tempdir = tempfile.mkdtemp(dir='.', prefix='tmp')
    for primer in primers:
        add_primer(primers=primer_dict, primer=primer)
    amplicons = list()
    if args.type == 'full' or args.type == 'no_primers':
        amplicons = create_amplicon(primers=primer_dict, type=args.type,
                                    path=tempdir)
    elif args.type == 'unique_amplicons':
        amplicons = create_unique_amplicon(primers=primer_dict,
                                           path=tempdir)
    else:
        sys.exit('Invalid amplicon type...')
    write_amplicon_to_bed(amplicons=amplicons, outfile=args.output)
    print(f"BED file written to {args.output}...")
    print(f"Removing temporary directory {tempdir}...")
    try:
        shutil.rmtree(tempdir)
    except:
        print(f"Unable to remove temporary directory {tempdir}")


if __name__ == '__main__':
    main()

