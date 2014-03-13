#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Chanjo - A clinical genomics coverage analysis tool

Usage:
  chanjo.py build <sql_db> [<reference>] [--dialect=STR] [--force]
  chanjo.py annotate <sql_db> <bam_file> [--cutoff=INT] [--extend-by=INT]
    [--dialect=STR] [--sample=STR] [--group=STR] [--contig-prepend=STR]
    [--threshold=INT] [--cores=INT] [--config=FILE] [--force]
  chanjo.py import <sql_db> [<annotate_output>] [--dialect=STR] [--json]
  chanjo.py read <bam_file> <contig> <start> <end> [--cutoff=INT]
    [--dialect=STR]
  chanjo.py peek <sql_db> <gene>... [--sample=STR] [--group=STR]
    [--dialect=STR]
  chanjo.py -h | --help
  chanjo.py -v | --version

Commands:
  build                Builds a new skeleton SQL interval store
  annotate             Annotate intervals in an existing database
  import               Import coverage annotations to an existing database
  read                 'Manually' read and calculate coverage for an interval
                       in a BAM-file
  peek                 Peek at annotated coverage metrics for supersets/genes

Arguments:
  <sql_db>             New or existing SQL database (path or URI)
  <bam_file>           Path to coverage source BAM-file e.g. './snap/align.bam'
  <reference>          Path to intervals reference dump e.g. CCDS database
  <gene>               List of HGNC gene symbols
  <annotate_output>    Output file from 'annotate' command (or stdin)
  <contig>             Contig (chromosome) ID e.g. chr2
  <start>              Starting interval position (1-based)
  <end>                Ending interval position (1-based)

Options:
  -h --help            Show this screen
  -v --version         Show version
  -c --cutoff=INT      Cutoff for completeness calculation [default: 10]
  -d --dialect=STR     Type of database: sqlite/mysql/BED [default: sqlite]
  -e --extend-by=INT   Extend intervals to include e.g. splice sites, +/- 2 bp
                       [default: 0]
  -s --sample=STR      Sample ID e.g. '0-0-0U' [default: default]
  -g --group=STR       Group ID e.g. to group samples in trios
                       [default: default]
  -t --threshold=INT   Base pair threshold for optimizing BAM-file reading
                       [default: 17000]
  -o --cores=INT       The maximum number of CPU cores to use
  -f --force           Overwrite existing files without warning
  -p --contig-prepend=STR   Prepend a string to each contig ID (USCS/Ensemble)
                       [default: ]
  -u --config=FILE     Path to YAML/JSON config file [default: global,local]
  -j --json            Import legacy JSON 'annotate' output file
"""

import json
import multiprocessing
import sys

from docopt import docopt
from clint.textui import puts, colored
from path import path

import chanjo
from chanjo.core import (annotate, build, import_data, import_json,
                         read_coverage, peek)


def open_or_stdin(file_path):
  """Opens a potential file or returns 'sys.stdin'

  Args:
    file_path (str): Path to potential file

  Returns:
    file: File handle for either the file or ``sys.stdin``
  """
  # Can't initialize ``path`` with ``None`` => point to non-existant path
  file_path = path(file_path or '__nonexistant')
  if file_path.isfile():
    return file_path.open('r')
  else:
    return sys.stdin


def main(args):
  # +------------------------------------------------------------------+
  # | Pre-process input arguments
  # +------------------------------------------------------------------+
  store_path = args['<sql_db>']
  coverage_path = args['<bam_file>']
  dialect = args['--dialect']
  sample_id = args['--sample']
  group_id = args['--group']
  FORCE = args['--force']
  cutoff = int(args['--cutoff'])

  # +------------------------------------------------------------------+
  # | Read an interval from a coverage source file (BAM)
  # +------------------------------------------------------------------+
  if args['read']:
    contig = args['<contig>']
    start = int(args['<start>'])
    end = int(args['<end>'])

    # Read and calculate coverage metrics
    _ = read_coverage(coverage_path, contig, start, end, cutoff)

  # +------------------------------------------------------------------+
  # | Peek at annotated coverage metrics for supersets/genes
  # +------------------------------------------------------------------+
  elif args['peek']:
    # Fetch the annotations from data store
    _ = peek(store_path, args['<gene>'], sample_id, dialect)

  # +------------------------------------------------------------------+
  # | Annotate intervals in an interval store
  # +------------------------------------------------------------------+
  elif args['annotate']:
    extension = int(args['--extend-by'])
    bp_threshold = int(args['--threshold'])
    prepend = args['--contig-prepend']

    # Use user defined number of cores or all that are available
    cores = int(args['--cores'] or multiprocessing.cpu_count())

    # Call main annotate function
    _ = annotate(sample_id, group_id, cutoff, coverage_path, store_path,
                 dialect, extension, prepend, bp_threshold)

  # +------------------------------------------------------------------+
  # | Import coverage annotations to an interval store
  # +------------------------------------------------------------------+
  elif args['import']:
    if args['--json']:
      # Import legacy JSON 'annotate' file with convertion to new format
      import_func = import_json
    else:
      import_func = import_data

    # Read from user defined file or piped standard in
    with open_or_stdin(args['<annotate_output>']) as input_stream:
      _ = import_func(store_path, input_stream, dialect)

  # +------------------------------------------------------------------+
  # | Builds a new skeleton SQL interval store
  # +------------------------------------------------------------------+
  elif args['build']:
    # Read CCDS reference file and import to a fresh SQL datastore
    _ = build(store_path, args['<reference>'], dialect, FORCE)

    puts(colored.blue('[chanjo]') + 'Built database: '
         + colored.green(store_path))


if __name__ == '__main__':
  # Parse docstring defined command line arguments
  args = docopt(__doc__, version='Chanjo v' + chanjo.__version__)

  # # Determine config file scope
  # scopes = args['--config'].split(',')

  # # If found, parse config files
  # options = rc.extend_args(args, __file__, defaults, scopes=scopes)

  main(args)
