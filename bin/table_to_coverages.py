#!/usr/bin/env python3
import argparse
import os
import sys
import polars as pl

def table_to_coverage_definitions(input_tsv, output_directory):
    tbl = pl.read_csv(input_tsv, separator = '\t')
    samples = tbl.columns[1:]
    for sample_id in samples:
        file = os.path.join(output_directory, f'{sample_id}.tsv')
        (
            tbl
            .select("genome_accession", sample_id)
            .filter(pl.col(sample_id) > 0)
            .write_csv(file, separator = "\t", include_header = False)
        )


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)

    parent_parser.add_argument('--profiles', help='Genome abundances per sample in tsv format', required=True)
    parent_parser.add_argument('--outdir', help='Output directory for coverages', required=True)

    args = parent_parser.parse_args()
    table_to_coverage_definitions(args.profiles, args.outdir)


