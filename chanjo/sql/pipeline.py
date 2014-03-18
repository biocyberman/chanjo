# -*- coding: utf-8 -*-

import itertools
import numpy as np

from .utils import load_ccds, filter_by_column, sort_by_columns
from .utils import process_superset


# +--------------------------------------------------------------------+
# | Importer pipeline
# +--------------------------------------------------------------------+
def import_from_ccds(db, ccds_path):
  # Build up the new schema
  db.setup()

  # Read in the entire CCDS database to a ``numpy.array``
  all_sets = load_ccds(ccds_path)

  # Filter by 'ccds_status' (we only care for 'Public' transcripts)
  public_sets = filter_by_column(all_sets, 'ccds_status', 'Public')

  # Sort the transcripts by gene (HGNC symbol) and chromosome to prepare
  # for grouping by gene.
  # We need sorting on 'chromosome' to make sure XY-genes are separated by
  # chromosome and don't get mixed up with each other.
  sorted_sets = sort_by_columns(public_sets, ('gene', 'chromosome'))

  # Group rows by gene/superset
  raw_supersets = itertools.groupby(sorted_sets, lambda x: x['gene'])

  # Extract and work on each of the groups of set/transcript rows
  processed_supersets = (process_superset(raw_superset)
                         for _, raw_superset in raw_supersets)

  # Store added intervals to avoid having to 'get_or_create' every time
  added_intervals = {}
  # Count how many superset we have processed and added
  count = 1
  # Work through all supersets etc.
  for superset_data, sets_and_intervals in processed_supersets:
    # Create a new superset and add it to the session
    new_superset = db.create('superset', *superset_data)
    db.add(new_superset)

    # Work through the sets in the superset
    for set_data, intervals in sets_and_intervals:
      # Create a new set
      new_set = db.create('set', *set_data)
      # Add the parent superset to the new set
      new_set.superset = new_superset
      # Add the new set to the session
      db.add(new_set)

      # Keep track of which interval is the first in the set
      first_interval = True
      # Work through each interval
      for interval_data in intervals:
        # Skip if we already processed+added this interval
        if interval_data[0] not in added_intervals:
          # Create new interval
          new_interval = db.create('interval', *interval_data)
          # Add the parent set to the interval => this also created the
          # intermediary 'interval__set' records
          new_interval.sets.append(new_set)

          # Mark each first interval as 'first'
          if first_interval:
            # First interval is first
            new_interval.first = True

          # Add new interval to the session
          db.add(new_interval)

          # Add new interval as 'added'
          added_intervals[interval_data[0]] = True

        # Turn off 'first_interval' flag
        first_interval = False

      # Last interval is last
      new_interval.last = True

    # Commit the session every 3000 supersets
    if count % 3000 == 0:
      db.commit()

    # Count each superset added
    count += 1

  # Commit the remaining supersets
  db.commit()

  return 0
