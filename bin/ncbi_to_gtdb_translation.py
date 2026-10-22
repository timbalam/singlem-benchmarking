#!/usr/bin/env python3
import argparse
import polars as pl

def read_ncbi_vs_gtdb_xlsx(path):
    d = (
        pl.read_excel(path, sheet_name = "Species",
                      columns = ["NCBI species", "List of R232 species"])
        .with_columns(
            pl.col("List of R232 species").str.split(', ')
        )
        .explode("List of R232 species")
        .with_columns(
            pl.col("List of R232 species").str.extract_groups(r"^(.+) ([^ ]+)%$")
        )
        .select(
            pl.col("NCBI species").alias("ncbi"),
            pl.col("List of R232 species").struct["1"].alias("r232"),
            pl.col("List of R232 species").struct["2"].cast(pl.Float64).alias("percent")
        )
        .with_columns(
            (pl.col("percent") == pl.col("percent").max().over(pl.col("ncbi")))
            .alias("which_max")
        )
        .filter(pl.col("which_max"))
        .drop("which_max")
    )
    assert d.shape[0] == d.select(pl.count("ncbi"))["ncbi"][0]
    return d

def write_ncbi_to_gtdb_mapping(mapping, output_file):
    mapping.write_csv(output_file, separator = '\t')

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)

    parent_parser.add_argument('--bac', help='NCBI vs GTDB r232 bacteria xlsx', required=True)
    parent_parser.add_argument('--ar', help='NCBI vs GTDB r232 archaea xlsx', required=True)
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
