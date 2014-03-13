# -*- coding: utf-8 -*-

import sys

import numpy as np


def get_chromosomes(prepend='chr', custom=[]):
  """Generate a list of human chromosome (contig) IDs.

  The `prepend` parameter allows for IDs matching both UCSC and NCBI
  standards (chr/no chr).

  Args:
    prepend (str, optional): What to prepend chromosome IDs (e.g. <chr>1)
    custom (list, optional): Additional chromosome/contig IDs to append.
      They will not be subject to the prepend option.

  Yields:
    str: The next chromosome ID in order, custom last
  """
  # Start with the autosomal chromosome IDs
  for num in range(1, 24):
    yield prepend + str(num)

  # Secondly yield the sex chromosomes
  for letter in ['X', 'Y']:
    yield prepend + letter

  # Lastly yield optional custom chromosome IDs without any prepend string
  for contig_id in custom:
    yield contig_id


def get_intervals(db, contig_id):
  """Returns all intervals for a given contig.

  Args:
    db (ElementAdapter): class:`chanjo.sql.ElementAdapter` class instance
    contig_id (str): Contig (chromosome) ID to match with

  Returns:
    list: List of intervals on the contig
  """
  interval = db.get('class', 'interval')
  return db.session.query(interval.start, interval.end, interval.id)\
           .filter_by(contig_id=contig_id).order_by(interval.start).all()


def group_intervals(intervals, threshold=1000, extension=0):
  """Groups and returns a list of intervals based on the threshold.

  Args:
    intervals (list): List of intervals
    threshold (int, optional): Approx. combined length per group
    extension (int, optional): Extend each interval +/- this number

  Yields:
    list of tuple: The next group of intervals
  """
  # Initialize stuff
  it = iter(intervals)
  interval = next(it)

  iStart = interval[0]
  iEnd = interval[1]
  interval_id = interval[2]

  # This is where we store grouped intervals + interval ID
  group = [(iStart - extension, iEnd + extension, interval_id)]

  for interval in it:

    start = interval[0]
    end = interval[1]
    interval_id = interval[2]

    # Optionally extend (widen) the segments
    start -= extension
    end += extension

    # Updated the current combined interval
    iEnd = max(iEnd, end)

    # If the current combined interval is big enough
    if (iEnd - iStart) > threshold:
      # Return currently grouped intervals
      yield group

      # Start a new combined interval
      group = [(start, end, interval_id)]

      # Reset the combined interval bounderies
      iStart, iEnd = start, end

    else:
      # Append to the current combined interval
      group.append((start, end, interval_id))

  # Return the last group
  yield group


def merge_intervals(intervals):
  """Returns the ends of a groups list of intervals

  Args:
    intervals (list): List of intervals

  Returns:
    tuple of int: The beginning and end of the combined interval
  """
  try:
    return intervals[0][0], intervals[-1][1]

  except IndexError:
    # The interval group didn't contain any intervals... why?
    return (0, 0)


def calculate_completeness(read_depths, threshold=10):
  """Calculates completeness across a number of positions.

  Note:
    The function catches the edge case where ``read_depths`` is an empty array
    which would lead to a ``ZeroDivisionError`` and returns 0% by default.

  Args:
    read_depths (array): :class:`numpy.array` of read depths for **each**
      of the positions
    threshold (int, optional): Cutoff to use for the filter

  Returns:
    float: The calculated completeness in percent
  """
  base_pair_count = len(read_depths)

  # Dodge rare division by zero error when `read_depths` is an empty array
  try:
    # Filter, then count bases with greater read depth than `threshold`
    # Divide by the total number of bases
    return len(read_depths[read_depths >= threshold]) / base_pair_count

  except ZeroDivisionError:
    # Without any bases to check, 0% pass the threshold
    return 0.


def assign_relative_positions(abs_start, abs_end, overall_start):
  """Return relative positions given the absolute positions and the overall
  starting position.

  Args:
    abs_start (int): Global start of the interval
    abs_end (int): Global end of the interval

  Returns:
    tuple of int: Relative start and end positions
  """
  rel_start = abs_start - overall_start
  rel_end = abs_end - overall_start

  return rel_start, rel_end


def calculate_values(read_depths, threshold):
  """Calculates mean coverage and completeness for a given read depths for
  a continous interval.

  Args:
    read_depths (array): :class:`numpy.array` of read depths for **each**
      of the positions
    threshold (int): Cutoff to use for the completeness filter

  Returns:
    tuple of float: Coverage and completeness for the interval represented by
      the read depth array.
  """
  return read_depths.mean(), calculate_completeness(read_depths, threshold)


def annotate_inverval_group(coverage_source, contig_id, interval_group,
                            cutoff):
  """Annotates and prints information about each segment in an interval group.

  Args:
    coverage_source (CoverageAdapter): Initialized
      :class:`chanjo.bam.CoverageAdapter` instance
    contig_id (str): Contig (chromosome) ID for all groups
    interval_group (list): List of intervals

  Returns:
    bool: ``True`` if successful, ``False`` otherwise.
  """
  # Pull through interval group to find out end of total interval
  listed_interval_group = list(interval_group)

  # Get the start and end of the group interval
  overall_start, overall_end = merge_intervals(listed_interval_group)

  # Get read depths for the whole (full) group interval
  read_depths = coverage_source.read(contig_id, overall_start, overall_end)

  # Loop through each of the intervals in the group
  for interval in listed_interval_group:
    # Convert to relative positions for the interval
    args = (interval[0], interval[1], overall_start)
    rel_start, rel_end = assign_relative_positions(*args)

    # Slice the overall read depth array with the relative coordinates
    read_depth_slice = read_depths[rel_start:rel_end + 1]

    # Calculate coverage and completeness for the interval
    coverage, completeness = calculate_values(read_depth_slice, cutoff)

    # Join the line components with a 'tab'-separator
    line = '\t'.join(map(str, [interval[2], coverage, completeness, '\n']))

    # Write to standard out (pipe)
    sys.stdout.write(line)

  return True
