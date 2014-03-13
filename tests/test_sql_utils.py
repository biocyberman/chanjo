# -*- coding: utf-8 -*-

import numpy as np

from .context import sql_utils


def test_read_supersets():
  # Dump the file into memory
  dump = np.recfromcsv('tests/fixtures/CCDS.mini.txt', delimiter='\t')

  # Filter by status code
  public_sets = dump[dump['ccds_status'] == b'Public']

  # Sort by HGNC symbol
  sorted_sets = public_sets[public_sets['gene'].argsort()]

  # Parse out the collected supersets
  supersets = list(sql_utils.read_supersets(sorted_sets))

  assert len(supersets) == 4
  assert len(supersets[0]) == 2
  assert len(supersets[1]) == 1
  assert supersets[0][0]['gene'] == b'RFPL2'
