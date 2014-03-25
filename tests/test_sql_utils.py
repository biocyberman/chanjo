# -*- coding: utf-8 -*-

import numpy as np

from .context import sql_utils

CCDS_PATH = 'tests/fixtures/CCDS.mini.txt'

RAW_SET = ['22', 'NC_000022.10', 'RFPL2', 10739, 'CCDS54521.1', 'Public',
           '-', 32586758, 32589173, '[32586758-32587338, 32588888-32589173]',
           'Identical']

RAW_SET_SEX = ['X', 'NC_000023.10', 'TRO', 7216, 'CCDS59529.1', 'Public', '+',
               54951423, 54957452, '[54951423-54951500]', 'Identical']


def test_load_ccds():
  # Dump the file into memory
  dump = sql_utils.load_ccds(CCDS_PATH)

  # Test some simple assertions
  assert len(dump) == 6
  assert len(dump[dump['chromosome'] == '22']) == 2

  # Test type conversions
  assert type(dump[0]['gene_id']) == np.int64


def test_filter_by_column():
  dump = sql_utils.load_ccds(CCDS_PATH)
  filtered_dump = sql_utils.filter_by_column(dump, 'gene', 'RFPL2')

  assert len(filtered_dump) == 2
  assert list(filtered_dump['chromosome']) == ['22', '22']


def test_sort_by_columns():
  dump = sql_utils.load_ccds(CCDS_PATH)

  # Test sorting by one column
  sorted_dump = sql_utils.sort_by_columns(dump, ('gene_id',))
  assert sorted_dump[0]['gene_id'] == 7216
  assert sorted_dump[-1]['gene_id'] == 728403


def test_parse_raw_intervals():
  one_interval = sql_utils.parse_raw_intervals('[801942-802433]')
  assert one_interval == [[801942, 802433]]

  multiple_intervals = sql_utils.parse_raw_intervals(
    "[9195451-9195936, 9196544-9196621, 9196750-9196861,"
    "9196963-9197108, 9197215-9197296, 9197991-9198013]")
  assert multiple_intervals == [[9195451, 9195936], [9196544, 9196621],
                                [9196750, 9196861], [9196963, 9197108],
                                [9197215, 9197296], [9197991, 9198013]]


def test_build_set():
  # Test random
  built_set = ('CCDS54521.1', '22', 32586758, 32589173, '-')
  assert sql_utils.build_set(RAW_SET) == built_set

  # Test sex chromosome set
  built_set = ('X-CCDS59529.1', 'X', 54951423, 54957452, '+')
  assert sql_utils.build_set(RAW_SET_SEX)


def test_build_interval():
  built_interval = sql_utils.build_interval((32586758, 32587338), RAW_SET)
  expected = ('22-32586758-32587338', '22', 32586758, 32587338, '-')
  assert built_interval == expected
