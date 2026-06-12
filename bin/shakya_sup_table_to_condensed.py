#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2025 Tim Lamberton
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

__author__ = "Tim Lamberton"
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

import polars as pl

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--supplementary-xlsx', help='Table S1 spreadsheet from Shakya et. al (2013)', required=True)
    parent_parser.add_argument('--sample', required=True)
    parent_parser.add_argument('--bac-tax', required=True)
    parent_parser.add_argument('--arc-tax', required=True)

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    taxonomy = (
        pl.concat([
            pl.read_csv(args.bac_tax, separator='\t', has_header=False, columns = 1, new_columns=['taxonomy']),
            pl.read_csv(args.arc_tax, separator='\t', has_header=False, columns = 1, new_columns=['taxonomy'])
        ])
        .unique()
        .with_columns(
            pl.col("taxonomy").str.extract("s__([^;]+)$").alias("Organism Name")
        )
    )
    pt = (
        pl.concat([
            pl.read_excel(args.supplementary_xlsx,
                          read_options = {"header_row": 2, "n_rows": 17,
                                          "use_columns": ["Organism Name", r"% of AB community genomes"]}),
            pl.read_excel(args.supplementary_xlsx,
                          read_options = {"header_row": 20, "n_rows": 49,
                                          "use_columns": ["Organism Name", r"% of AB community genomes"]})
        ])
        .rename({r"% of AB community genomes": "coverage"})
        .join(taxonomy, on = "Organism Name", how = "left",
              validate = "1:1")
        .sort(pl.col("taxonomy"))
        .filter(pl.col("coverage") > 0)
    )
    import pdb; pdb.set_trace()
    output = (
        pt.select([
            pl.lit(args.sample).alias("sample"),
            pl.col("coverage"),
            pl.col("taxonomy")
        ])
        .write_csv(include_header=True, separator='\t')
    )
    print(output, end="")

    logging.info("Done")

