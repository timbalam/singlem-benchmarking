#!/usr/bin/env python3
import argparse
import polars as pl


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)

    parent_parser.add_argument('--profile', help='NCBI vs GTDB r232 bacteria xlsx', required=True)
    parent_parser.add_argument('--output-file', help='Output mapping table', required=True)

    args = parent_parser.parse_args()
    mapping = (
        pl.concat([
            read_ncbi_vs_gtdb_xlsx(args.bac)
            .with_columns(domain = pl.lit("bacteria")),
            read_ncbi_vs_gtdb_xlsx(args.ar)
            .with_columns(domain = pl.lit("archaea"))
        ])
        .sort('percent')
    )
    write_ncbi_to_gtdb_mapping(mapping, args.output_file)
