#! /usr/bin/env python

import argparse
import os
import sys
import logging
import extern
import tempfile
import polars as pl
import shutil

def simulate_art(*, read_length, coverages, genomes, taxonomy, output1,
                 output2, art_bin, threads):

    df = (
        coverages
        .join(taxonomy, on = 'taxonomy', how = 'inner')
        .join(genomes, on = 'otu', how = 'inner')
        .with_columns(
            cmd = pl.lit(shutil.which(art_bin))
                + ' -ss HSXt -i '
                + pl.col("path")
                + ' -p -l '
                + pl.format("{}", pl.lit(read_length))
                + ' -f '
                + pl.format("{}", pl.col("coverage"))
                + ' -m 400 -s 10 -o simulated_reads/'
                + pl.format("{}", pl.int_range(pl.len()))
                + '. &> /dev/null'
        )
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        os.makedirs('simulated_reads')
        sim_commands = df["cmd"]
        logging.info(f"Simulating {len(sim_commands)} genomes ..")
        extern.run_many(sim_commands, num_threads=threads, progress_stream=sys.stderr)

        logging.info("Concatenating simulated reads and compressing ..")
        extern.run("cat simulated_reads/*1.fq |sed 's=/= =' |pigz -p {} >{}".format(threads, output1))
        extern.run("cat simulated_reads/*2.fq |sed 's=/= =' |pigz -p {} >{}".format(threads, output2))

def read_coverage_file(path, sample):
    return (
        pl
        .read_csv(path, separator = '\t', columns = ['sample', 'taxonomy', 'coverage'])
        .filter(pl.col("sample") == sample)
    )

def read_taxonomy_file(path):
    return pl.read_csv(path, separator = '\t', has_header = False,
                       new_columns = ["otu", "taxonomy"])

def read_genomes_list(path):
    dir = os.path.dirname(path)
    make_absolute = lambda path: os.path.normpath(os.path.join(dir, path))
    return (
        pl.read_csv(path, separator = '\t', has_header = False,
                    new_columns = ["path", "otu"])
        .with_columns(
            pl.col("path").map_elements(make_absolute)
        )
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument('--coverage-file', required=True, help='Path to file with taxonomic + coverage information (tsv with sample/taxonomy/coverage)')
    parser.add_argument('--genome-list', required=True, help='Path to file with OTU names and genome file paths (tsv with otu/path)')
    parser.add_argument('--taxonomy', required=True, help='Path to file with taxonomies for otus in genomes-list (tsv with otu/taxonomy)')
    parser.add_argument('-1', '--read1', required=True, help='Path to output fq.gz file')
    parser.add_argument('-2', '--read2', required=True, help='Path to output fq.gz file')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('--art-bin', required=True, help='Path to ART binary (e.g. art_illumina)')
    parser.add_argument('--read-length', type=int, default = 150, help='Simulated read length')
    parser.add_argument('--sample', help='Name of sample (in coverage file)')
    
    args = parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    
    simulate_art(
        read_length = args.read_length,
        output1 = os.path.abspath(args.read1),
        output2 = os.path.abspath(args.read2),
        coverages = read_coverage_file(args.coverage_file, args.sample),
        genomes = read_genomes_list(args.genome_list),
        taxonomy = read_taxonomy_file(args.taxonomy),
        threads = args.threads,
        art_bin = args.art_bin
    )
    