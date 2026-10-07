#!/usr/bin/env python3
import argparse
import os
import sys
from cami_opal.utils import load_data
sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

def biobox_to_condensed_profiles(input_biobox_path, output_directory):
    sample_list = load_data.open_profile_from_tsv(input_biobox_path, False)
    for sample_id, header, profile in sample_list:
        file = os.path.join(output_directory, f'{sample_id}.tsv')

        with open(file, 'w') as f:
            f.write("\t".join(["sample", "coverage", "taxonomy"]) + "\n")
            for prediction in profile:
                if prediction.rank == "strain" and prediction.percentage > 0:
                    f.write("\t".join([sample_id, str(prediction.percentage), prediction.taxid]) + "\n")


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)

    parent_parser.add_argument('--profiles', help='CAMI profile to extract taxids from', required=True)
    parent_parser.add_argument('--outdir', help='Output directory for condensed profiles', required=True)

    args = parent_parser.parse_args()
    biobox_to_condensed_profiles(args.profiles, args.outdir)


