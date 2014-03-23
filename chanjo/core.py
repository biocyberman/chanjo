# -*- coding: utf-8 -*-
"""
chanjo.core
~~~~~~~~~~~

Contains the individual pipelines that are invoked from the command line
script.
"""

import errno
import json
import sys
import csv

from clint.textui import puts, colored
from path import path
from sqlalchemy.exc import IntegrityError

from .bam import CoverageAdapter
from .sql.core import ElementAdapter
from .sql.models import Base, Sample, SupersetData
from .sql.pipeline import import_from_ccds
from .utils import get_chromosomes, get_intervals, group_intervals
from .utils import annotate_inverval_group, calculate_values


def annotate(sample_id, group_id, cutoff, bam_path, sql_path, dialect,
             extension, prepend, bp_threshold):
  """Automates a pipeline for annotating all intervals from a SQL database.
  Writes both metadata (header) and tabular data for each interval
  with calculated coverage and completeness to standard out.

  Args:
    sample_id (str): Unique Id for a given sample
    group_id (str): Id for a given group of samples (e.g. family/trio)
    cutoff (int): Threshold to use for completeness calculation
    bam_path (str): Path to BAM-file
    sql_path (str): Path to existing ('built') SQL database
    dialect (str): SQL database dialect (+ optionally Python adapter to use)
    extension (int): Number of bases to extend each interval with (+/-)
    prepend (str): Renames each contig by prepending this string
    bp_threshold (int): Optimization number for reading BAM-file in chunks

  Returns:
    bool: Exit code (success or fail)
  """
  # Write metadata to output header
  header = {
    'sample_id': sample_id,
    'group_id': group_id,
    'cutoff': cutoff,
    'coverage_source': path(bam_path).abspath(),
    'element_source': path(sql_path).abspath(),
    'extension': extension
  }
  line = '#{}\n'.format(json.dumps(header))

  # Print to standard out (pipe)
  sys.stdout.write(line)

  # Connect to ElementAdapter and CoverageAdapter
  db = ElementAdapter(sql_path, dialect=dialect)
  coverage = CoverageAdapter(bam_path)

  # Generate list of contig Ids
  contig_ids = list(get_chromosomes(prepend=prepend))

  # 'Exchange' contig Ids for list of contig intervals (1-based positions)
  contigs = (get_intervals(db, contig_id) for contig_id in contig_ids)

  # 'Exchange' contig intervals for list of grouped contig intervals
  # As for positions, 1/0-based doesn't matter, only the lenght is considered
  grouped_intervals = (group_intervals(intervals, bp_threshold, extension)
                       for intervals in contigs)

  # Process a group of intervals
  # 1. Read from coverage source for each full group of intervals
  # 2. For each interval in each group: calculate coverage and completeness
  # 3. Print to stdout for each processed interval
  for contig_id, contig_group in zip(contig_ids, grouped_intervals):
    for group in contig_group:
      annotate_inverval_group(coverage, contig_id, group, cutoff)

  return 0


def import_data(sql_path, input_stream, dialect):
  """Import output from 'annotate' sub-command.

  Args:
    sql_path (str): Path to the database or MySQL connection string
    input_stream (file): File handle (file or stdin)
    dialect (str): SQL database type (sqlite or mysql+<connector>)

  Returns:
    bool: True if successful
  """
  db = ElementAdapter(sql_path, dialect=dialect)

  # Load the header information as JSON (removing the '#' character)
  metadata = json.loads(input_stream.readline()[1:])

  # Extract info used later on
  sample_id = metadata['sample_id']
  group_id = metadata['group_id']

  # Add a Sample entry with metadata
  db.add(db.create('sample', **metadata))

  # Loop through the rest of the lines which should be tab–separated
  for line in csv.reader(input_stream, delimiter='\t'):

    # Create a new intervals data entry
    # Args: parent_id, coverage, completeness, ...
    args = (line[0], sample_id, group_id, float(line[1]), float(line[2]))
    db.add(db.create('interval_data', *args))

  # Commit updates after loading all intervals
  db.commit()

  # Extend annotations to sets and supersets
  return extend_annotations(db, sample_id, group_id)


def import_json(sql_path, input_stream, dialect):
  """Legacy importer for the old JSON output format of the 'annotate' command.
  """
  db = ElementAdapter(sql_path, dialect=dialect)

  dump = json.load(input_stream)

  annotations = dump['annotations']
  sample_id = dump['sample_id']
  group_id = dump['group_id']

  if dump['splice']:
    extension = 2
  else:
    extension = 0

  # Add a Sample entry with metadata
  sample = db.create(
    'sample',
    sample_id=sample_id,
    group_id=str(group_id),
    cutoff=dump['cutoff'],
    extension=extension,
    coverage_source=dump['source'],
    element_source=None
  )
  db.add(sample)

  # For each of the annotations (intervals)
  for annotation in annotations:
    interval_data = db.create(
      'interval_data',
      parent_id=_convert_old_interval_id(annotation[2]),
      coverage=annotation[0],
      completeness=annotation[1],
      sample_id=None,
      group_id=group_id
    )
    interval_data.sample = sample
    db.add(interval_data)

  # Commit intervals before proceeding
  # We do this in part to leverage subsequent SQL queries.
  db.commit()

  # Extend annotations to sets and supersets
  return extend_annotations(db, sample_id, group_id)


def _convert_old_interval_id(old_id):
  """Private function for converting a 0:0-based exon Id to the new 1:1-based
  interval Id.

  Args:
    old_id (str): Old exon Id, '0:0-based'

  Returns:
    str: New interval Id, '1:1-based'
  """
  # Split into parts (contig, start, end)
  parts = old_id.split('-')

  # Recombine but with converted coordinates from 0:0 to 1:1
  return '-'.join([parts[0], str(int(parts[1]) + 1), str(int(parts[2]) + 1)])


def extend_annotations(db, sample_id, group_id):
  """Extend interval annotations to sets and supersets by calculating the mean
  of related elements.

  Args:
    db (CoverageAdapter): Instance of class:`chanjo.sql.core.CoverageAdapter`
    sample_id (str): ID of sample to extend
    group_id (str): Group ID of sample to extend

  Returns:
    bool: ``True`` if successful, ``False`` otherwise.
  """
  # Extend interval annotations to sets
  # We commit so the SQL query in the next step works
  db.add([db.create(
    'set_data',
    parent_id=raw_set[0],
    sample_id=sample_id,
    group_id=group_id,
    coverage=raw_set[1],
    completeness=raw_set[2]
  ) for raw_set in db.set_stats(sample_id)]).commit()

  # Extend annotations to genes
  db.add([db.create(
    'superset_data',
    parent_id=raw_superset[0],
    sample_id=sample_id,
    group_id=group_id,
    coverage=raw_superset[1],
    completeness=raw_superset[2]
  ) for raw_superset in db.superset_stats(sample_id)]).commit()

  return True


def build(sql_path, ccds_path, dialect, force=False):
  """Builds up a new database skeleton or wipes an existing one. As reference
  it uses the CCDS database dump of transcript annotations.

  Args:
    sql_path (str): Path to the datastore to setup
    ccds_path (str): Path to the reference CCDS database
    dialect (str): SQL dialect
    force (bool, optional): Option to wipe an existing before rebuilding
        without warnings. Defaults to ``False``.

  Returns:
    bool: ``True`` if successful, ``False`` otherwise.
  """
  # Connect to the database
  db = ElementAdapter(sql_path, dialect=dialect)

  # Check if the database already exists (expect 'mysql' to exist)
  # 'dialect' is in the form of '<db_type>+<connector>'
  if dialect.startswith('mysql') or path(sql_path).exists():
    if force:
      # Wipe the database clean with a warning
      puts(colored.blue('[chanjo] ') + colored.red('Wiping database...'))
      db.tare_down()
    elif dialect == 'sqlite':
      # Prevent from wiping existing database to easily
      raise OSError(errno.EEXIST, sql_path)

  # Set up new tables with a notice
  puts(colored.blue('[chanjo]') + ' ' + colored.green('Building database...'))
  db.setup()

  try:
    # Start the import
    _ = import_from_ccds(db, ccds_path)
  except IntegrityError:
    puts("Database was already built! Please use '--force' to overide "
         "this error")
    return False

  return True


def read_coverage(bam_path, contig_id, start, end, cutoff):
  """Takes a genomic interval and calculates coverage and completeness metrics.
  Writes results to ``sys.stdout``; coverage and completeness.

  Args:
    bam_path (str): Path to the coverage source BAM-file
    contig_id (str): Contig/chromosome ID
    start (int): Beginning of interval (1-based)
    end (int): End of interval (1-based)
    cutoff (int): Lower threshold for completeness calculation

  Returns:
    bool: ``True`` if successful, ``False`` otherwise.
  """
  try:
    # Connect to the coverage source (BAM-file)
    coverage_source = CoverageAdapter(bam_path)
  except OSError:
    sys.exit("The file doesn't exist: {}".format(bam_path))

  # Write header for tabular output
  sys.stdout.write('#coverage\t#completeness\n')

  # Read read depths from the coverage source
  read_depths = coverage_source.read(contig_id, start, end)

  # Calculate mean coverage and completeness for the interval
  coverage, completeness = calculate_values(read_depths, cutoff)

  # Write the output to stdout
  sys.stdout.write('{0}\t{1}\n'.format(coverage, completeness))

  return True


def peek(sql_path, superset_ids, sample_id=None, dialect='sqlite'):
  """Get a peak at some of the annotations (supersets) in a database. Writes
  output to ``sys.stdout`` (one line per sample and superset).

  Args:
    sql_path (str): Path to the SQL database
    superset_ids (list): List of superset IDs
    sample_id (str, optional): Limit results to a single sample

  Returns:
    bool: ``True`` if successful, ``False`` otherwise.
  """
  # Connect to intervals store
  db = ElementAdapter(sql_path, dialect=dialect)

  # Write header for tabular output
  sys.stdout.write('#superset\tsample\tcoverage\tcompleteness\n')

  # Collect the requested data
  output = {}

  # Build base query
  sd = SupersetData
  query = db.query(sd.parent_id, sd.sample_id, sd.coverage, sd.completeness)\
            .filter(sd.parent_id.in_(superset_ids))

  # Optionally filter the query by sample ID
  if sample_id:
    query = query.filter(SupersetData.sample_id == sample_id)

  # Loop over each of the query results
  for result in query.all():

    # Write output to stdout
    line = '\t'.join([result[0], result[1], str(result[2]), str(result[3])])
    sys.stdout.write(line + '\n')

  return True
