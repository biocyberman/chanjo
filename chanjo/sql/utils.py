# -*- coding: utf-8 -*-


def read_supersets(sorted_sets, superset_column='gene'):
  """Yields collections of sets for each superset in the input. Expects the
  input to be sorted by the 'gene' column.
  """
  # Initialize with the first superset
  last_superset = sorted_sets[0][superset_column]
  sets = []

  # Record is a CCDS row, representing a set
  for record in sorted_sets:

    # Check if the record belongs to the same superset as last
    if record[superset_column] == last_superset:
      # Add the set to the collection of sets
      sets.append(record)

    else:
      # Yield the sets now representing a whole superset
      yield sets
      # Reset the sets with the currect (new) set
      sets = [record]
      # Update last superset
      last_superset = record[superset_column]

  # Also yield the last collection of sets for the last superset
  yield sets


def parse_raw_intervals(str_list):
  """Takes the formatted string of interval coordinates from the CCDS row
  and turns it into a more managable list of lists with (start, end)
  coordinates for each interval.

  :param str str_list: A csv string of (start,end) pairs, wrapped in '[]'
  :returns: A list of lists with the start, end pairs
  :rtype: int, int
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


def assemble_interval(raw_set, raw_interval):
  contig_id = raw_set['chromosome']

  return (
    '{}-{}-{}'.format(contig_id, *raw_interval),
    contig_id,
    raw_interval[0],
    raw_interval[1],
    raw_set['cds_strand']
  )


def assemble_set(raw_set):
  contig_id = raw_set['chromosome']
  set_id = raw_set['ccds_id']

  # 20 supersets are present on both X/Y (this trickles down to sets)
  if raw_set['chromosome'] == 'X' or raw_set['chromosome'] == 'Y':
    set_id = '{}-{}'.format(contig_id, set_id)

  set_data = (
    set_id,
    contig_id,
    raw_set['cds_from'],
    raw_set['cds_to'],
    raw_set['cds_strand'],
    raw_set['gene']
  )

  interval_coordinates = parse_raw_intervals(raw_set['cds_locations'])

  intervals = [assemble_interval(raw_set, raw_interval)
               for raw_interval in interval_coordinates]

  interval_sets = [('{}-{}-{}'.format(contig_id, *raw_interval), set_id)
                   for raw_interval in interval_coordinates]

  return set_data, intervals, interval_sets


def assemble_superset(raw_superset):
  first_set = raw_superset[0]

  contig_id = first_set['chromosome']
  superset_id = first_set['gene']

  # 20 supersets are present on both X/Y
  if first_set['chromosome'] == 'X' or first_set['chromosome'] == 'Y':
    superset_id = '{}-{}'.format(contig_id, superset_id)

  superset = (
    superset_id,
    contig_id,
    min([raw_set['cds_from'] for raw_set in raw_superset]),
    max([raw_set['cds_to'] for raw_set in raw_superset]),
    first_set['cds_strand'],
    first_set['gene_id']
  )

  sets_data = [assemble_set(raw_set) for raw_set in raw_superset]

  return superset, sets_data
