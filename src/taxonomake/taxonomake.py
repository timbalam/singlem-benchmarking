#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2021 Ben Woodcroft, 2025 Tim Lamberton
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
import os
import tempfile
from taxonomake.modules.community_description import process_community_description

def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parser.add_argument('--quiet', help='only output errors', action="store_true")
    parser.add_argument('-d', '--directory', help='Directory for intermediate outputs (defaults to a temporary directory)', dest='directory',
                        default=None)
    parser.add_argument(
        '--snakemake-args',
        help='Additional arguments to supplied to snakemake in the form of a single string '
             'e.g. "--print-compilation True". \n '
             'NOTE: Most commands in snakemake -h are valid but some commands may clash with commands \n '
             'taxonomake directly supplies to snakemake. Please make sure your additional commands don\'t clash.',
        default='',
    )
    parser.add_argument('configfile', default = "community.toml")

    args = parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    prefix = args.directory
    if prefix is not None:
        os.makedirs(prefix, exist_ok = True)
        def cleanup():
            pass
    else:
        td = tempfile.TemporaryDirectory()
        prefix = td.name
        def cleanup():
            td.cleanup()
    try:
        process_community_description(args.configfile, prefix = prefix,
                                      snakemake_args = args.snakemake_args)
    finally:
        cleanup()

if __name__ == '__main__':
    main()

