# -*- coding: utf-8 -*-

import numpy as np

from .context import sql_utils


def test_read_supersets():
  # Dump the file into memory
  ccds_path = 'tests/fixtures/CCDS.mini.txt'
  columns = [('chromosome', '<U12'), ('nc_accession', '<U12'),
             ('gene', '<U24'), ('gene_id', int), ('ccds_id', '<U12'),
             ('ccds_status', '<U12'), ('cds_strand', '<U12'),
             ('cds_from', int), ('cds_to', int), ('cds_locations', '<U10000'),
             ('match_type', '<U12')]
  dump = np.genfromtxt(ccds_path, delimiter='\t', dtype=columns)

  # Filter by status code
  public_sets = dump[dump['ccds_status'] == 'Public']

  # Sort by HGNC symbol
  sorted_sets = public_sets[public_sets['gene'].argsort()]

  # Parse out the collected supersets
  supersets = list(sql_utils.read_supersets(sorted_sets))

  assert len(supersets) == 4
  assert len(supersets[0]) == 2
  assert len(supersets[1]) == 1
  assert supersets[0][0]['gene'] == 'RFPL2'
