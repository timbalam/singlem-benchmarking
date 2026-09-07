#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2021 Ben Woodcroft
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2021"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
import logging
import sys
import os

import pandas as pd
import polars as pl

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--genome-to-id', help='Mapping from genome to ID', required=True)
    parent_parser.add_argument('--genome-lengths-directory', help='Directory of genome lengths', required=True)
    parent_parser.add_argument('--read-mapping-counts', help='Read mapping counts', required=True)

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Get the genome lengths of each CAMI genome
    genome_to_length = {}
    for tsv in os.listdir(args.genome_lengths_directory):
        if 'tsv' not in tsv:
            continue
        gdf = pl.read_csv(os.path.join(args.genome_lengths_directory, tsv), separator='\t', has_header=False)
        cami_genome_name = os.path.basename(tsv).replace('.fasta.tsv','')
        genome_length = gdf['column_2'].sum()
        genome_to_length[cami_genome_name] = genome_length
    logging.info("Read in %d genome lengths" % len(genome_to_length))

    # Read genome to Id mapping
    # Otu520.0        /net/sgi/cami/data/CAMI2/nwillassen_data/marine/simulation_short_read/genomes/Filomicrobium_sp_W1_genomic.fasta
    # Otu1947 /net/sgi/cami/data/CAMI2/nwillassen_data/marine/simulation_short_read/genomes/Shewanella_sp_MR-4_genomic.fasta
    # Otu806.0        /net/sgi/cami/data/CAMI2/nwillassen_data/marine/simulation_short_read/genomes/Shewanella_violacea_DSS12_genomic.fasta
    genome_to_id_df = pl.read_csv(args.genome_to_id, separator='\t', has_header=False)
    # Want OTU to basename minus .fasta
    genome_to_id = {}
    for row in genome_to_id_df.rows():
        genome_to_id[row[0]] = os.path.basename(row[1]).replace('.fasta','')

    mapping_counts = (
        pl.read_csv(args.read_mapping_counts, separator=' ', has_header=True, new_columns = ['read_count', 'otu'])
    )
    mapping_counts = (
        mapping_counts
        .with_columns(
            cami_name = pl.col('otu').replace(genome_to_id).alias('cami_name')
        )
        .with_columns(
            genome_length = pl.col('cami_name').replace_strict(genome_to_length)
        )
        .with_columns(
            coverage = pl.col('read_count') * 150 / pl.col('genome_length')
        )
    )

    (
        mapping_counts
        .select(['otu', 'coverage'])
        .write_csv(sys.stdout, separator = '\t', include_header = False)
    )


