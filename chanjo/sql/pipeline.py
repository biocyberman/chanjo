# -*- coding: utf-8 -*-

import numpy as np

from .core import ElementAdapter
from .utils import read_supersets, assemble_superset


# +--------------------------------------------------------------------+
# | Import pipeline
# +--------------------------------------------------------------------+
def import_from_ccds(db, ccds_path):
  # Dump the file into memory
  all_sets = np.recfromcsv(ccds_path, delimiter='\t')

  # Filter out all non-Y records
  non_y = all_sets[all_sets['chromosome'] != b'Y']
  all_y = all_sets[all_sets['chromosome'] == b'Y']

  for full_set in [non_y, all_y]:

    # Filter by status code
    public_sets = full_set[full_set['ccds_status'] == b'Public']

    # Sort by HGNC symbol
    sorted_sets = public_sets[public_sets['gene'].argsort()]

    # Yield each gene
    raw_supersets = read_supersets(sorted_sets)

    # Process each superset
    supersets = map(assemble_superset, raw_supersets)

    count = 0
    added_intervals = {}
    added_sets = {}
    for superset_data, sets in supersets:
      new_superset = db.create('superset', *superset_data)

      for set_data, intervals, interval_sets in sets:
        if set_data[0] not in added_sets:
          added_sets[set_data[0]] = True
          new_set = db.create('set', *set_data)
        else:
          print(set_data)

        first_interval = True
        for interval_data in intervals:
          if interval_data not in added_intervals:
            added_intervals[interval_data] = True
            new_interval = db.create('interval', *interval_data)
            new_set.intervals.append(new_interval)

            if first_interval:
              # First interval is first
              new_interval.first = True
              # ... and no other intervals
              first_interval = False

            db.add(new_interval)

        # Last interval is last
        new_interval.last = True

        db.add(new_set)

      db.add(new_superset)

      if count == 5000:
        db.commit()
        count = 0

      else:
        count += 1

    # Commit remaining entries
    db.commit()

    return 0
