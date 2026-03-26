#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2026 Tim Lamberton
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
__copyright__ = "Copyright 2026"
__credits__ = ["Tim Lamberton"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
import logging
#import sys
import os
import json
from random import shuffle

#sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input-archive-otu-table', required=True)
    parent_parser.add_argument('--fold', required=True)
    parent_parser.add_argument('--output-mask-dir', required=True)

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    output_prefix = f"{args.output_mask_dir}/{os.path.splitext(os.path.basename(args.input_archive_otu_table))[0]}.mask"

    with open(args.input_archive_otu_table) as io:
        j = json.load(io)
    seq_field = j['fields'].index('sequence')

    shuffle_otus = [o[seq_field] for o in j['otus']]
    shuffle(shuffle_otus)
    
    n_fold = int(args.fold)
    n_split = round(len(shuffle_otus) / n_fold)
    for i in range(n_fold):
        start = i * n_split
        stop = min((i + 1) * n_split, len(shuffle_otus))
        split_otus = shuffle_otus[start:stop]

        with open(f"{output_prefix}{i}.txt", "w") as f:
            for seq in split_otus:
                f.write(f"{seq}\n")
        
    logging.info("Done")
