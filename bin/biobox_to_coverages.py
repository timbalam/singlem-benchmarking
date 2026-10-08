#!/usr/bin/env python3
import argparse
import os
import sys
from cami_opal.utils import load_data
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

def biobox_to_coverage_definitions(input_biobox_path, output_directory):
    #https://github.com/CAMI-challenge/OPAL/blob/master/cami_opal/utils/load_data.py
    sample_list = load_data.open_profile_from_tsv(input_biobox_path, False)
    for sample_id, header, profile in sample_list:
        file = os.path.join(output_directory, f'{sample_id}.tsv')

        with open(file, 'w') as f:
            for prediction in profile:
                if prediction.rank == "strain" and prediction.percentage > 0:
                    f.write("\t".join([prediction.taxid, str(prediction.percentage)]) + "\n")


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)

    parent_parser.add_argument('--profiles', help='CAMI profile to extract taxids from', required=True)
    parent_parser.add_argument('--outdir', help='Output directory for condensed profiles', required=True)

    args = parent_parser.parse_args()
    biobox_to_coverage_definitions(args.profiles, args.outdir)


