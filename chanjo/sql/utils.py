# -*- coding: utf-8 -*-

import numpy as np


def load_ccds(path):
  """Read a CCDS database dump into a numpy 2D array with correct convertion to
  unicode string and integers.

  Args:
    path (str): Path to the database dump

  Returns:
    numpy.array: 2D array with defined column names
  """
  # Set up column names and types (otherwise we get bytes in Python 3)
  columns = [('chromosome', '<U12'), ('nc_accession', '<U12'),
             ('gene', '<U24'), ('gene_id', int), ('ccds_id', '<U12'),
             ('ccds_status', '<U12'), ('cds_strand', '<U12'),
             ('cds_from', int), ('cds_to', int), ('cds_locations', '<U7700'),
             ('match_type', '<U12')]

  # Dump the file into memory
  return np.genfromtxt(path, delimiter='\t', dtype=columns)


def filter_by_column(array, column_id, value):
  """Filter ``numpy.array`` by a value in a column.

  Args:
    array (numpy.array): Array to filter
    column_id (str): Id for column to filer on
    value (any): Usually a string to match against

  Returns:
    numpy.array: Subset of rows with matching column values
  """
  return array[array[column_id] == value]


def sort_by_columns(array, column_ids):
  """Sort a ``numpy.array`` by multiple column Ids.

  Args:
    array (numpy.array): Array to sort
    column_ids (tuple): Column Ids to sort on (note that order matters)

  Returns:
    numpy.array: Sorted array
  """
  return array[np.argsort(array, order=column_ids)]


def build_superset(raw_superset):
  """Extract 'superset data' from a subset of ``numpy.array`` records
  (CCDS database rows) with common superset Id.

  Args:
    raw_superset (numpy.array): Set of records with common superset Id

  Returns:
    tuple: All data needed to create a new overall 'superset' record
  """
  # Extract information that *should* be the same for each of the child sets
  any_set = raw_superset[0]
  contig_id = any_set['chromosome']
  superset_id = any_set['gene']

  # 20 supersets are present on both X/Y
  if contig_id in ('X', 'Y'):
    superset_id = '{}-{}'.format(contig_id, superset_id)

  return (
    superset_id,
    contig_id,
    int(min([raw_set['cds_from'] for raw_set in raw_superset])),
    int(max([raw_set['cds_to'] for raw_set in raw_superset])),
    any_set['cds_strand'],
    int(any_set['gene_id'])
  )


def build_set(raw_set):
  """Extract 'set data' from a single ``numpy.array`` record (CCDS row).

  Args:
    raw_set (tuple): Row from CCDS database (numpy.void also works)

  Returns:
    tuple: All data needed to create a new 'set' record
  """
  # Extract information that *should* be the same for each of the child sets
  contig_id = raw_set[0]
  set_id = raw_set[4]

  # 20 supersets are present on both X/Y
  if contig_id in ('X', 'Y'):
    set_id = '{}-{}'.format(contig_id, set_id)

  return (
    set_id,
    contig_id,
    int(raw_set[7]),  # From
    int(raw_set[8]),  # To
    raw_set[6]        # Strand
  )


def build_interval(raw_interval, raw_set):
  """Extract 'interval data' for a single interval (exon).

  Args:
    raw_interval (tuple): Pair of start and end coordinates for the interval
    raw_set (tuple): Relevant row from CCDS database (numpy.void also works)

  Returns:
    tuple: All data needed to create a new 'interval' record
  """
  contig_id = raw_set[0]

  return (
    '{}-{}-{}'.format(contig_id, *raw_interval),
    contig_id,
    int(raw_interval[0]),
    int(raw_interval[1]),
    raw_set[6]        # Strand
  )


def process_superset(raw_superset):
  """Process a group of rows from the CCDS database with common superset Id.

  Args:
    raw_superset (generator): ``itertools.groupby`` generator

  Returns:
    tuple: superset data and generator for yielding set and interval data
  """
  # Unfold and store the generator so we can use the result more than once
  raw_superset = list(raw_superset)
  # Extract superset data
  superset_data = build_superset(raw_superset)
  # 'Recursively' mine each set/row for set and interval data
  sets_and_intervals = (process_set(raw_set) for raw_set in raw_superset)

  return superset_data, sets_and_intervals


def process_set(raw_set):
  """Process a row from the CCDS database.

  Args:
    raw_superset (tuple): Row from the CCDS database with unique set Id

  Returns:
    tuple: set data and generator for yielding interval data
  """
  # Extract data for the set
  set_data = build_set(raw_set)
  # Parse the 'cds_location' string and mine each pair of (start, end)
  # coordinates for interval data
  intervals = (build_interval(raw_interval, raw_set)
               for raw_interval in parse_raw_intervals(raw_set[9]))

  return set_data, intervals


def parse_raw_intervals(str_list):
  """Takes the formatted string of interval coordinates from the CCDS row
  and turns it into a more managable list of lists with (start, end)
  coordinates for each interval.

  Args:
    str_list (str): A csv string of (start,end) pairs, wrapped in '[]'

  Returns:
    list: 2D list with the start ``int``, end ``int`` pairs
  """
  # Remove the "[]"
  csv_intervals = str_list[1:-1].replace(' ', '')

  # 1. Split first into exons coordinates
  # 2. Split into start, end and parse int
  intervals = [[int(pos) for pos in item.split('-')]
               for item in csv_intervals.split(',')]

  # 3. Correct coords to 0,0-based Pythonic standard
  for interval in intervals:
    interval[0] -= 1
    interval[1] -= 1

  return intervals
