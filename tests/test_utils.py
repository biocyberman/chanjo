# -*- coding: utf-8 -*-

import sys

import numpy as np

from .context import utils, ElementAdapter


def test_id_generator():
  # Test getting different length random Ids
  ten_char_id = utils.id_generator(10)
  assert len(ten_char_id) == 10

  # Test with custom input chars
  custom_id = utils.id_generator(chars='ATCG')
  assert 'U' not in custom_id


def test_open_or_stdin():
  # Test with existing file (expect file handle)
  existing_file = 'tests/fixtures/CCDS.mini.txt'
  f = utils.open_or_stdin(existing_file)
  with open(existing_file, 'r') as handle:
    assert f.name == handle.name
  f.close()

  # Test with non-existent file (expect stdin)
  nonexisting_file = 'tests/fixtures/i_dont_exist.txt'
  stdin = utils.open_or_stdin(nonexisting_file)
  assert stdin == sys.stdin


def test_get_chromosomes():
  # Test with empty prepend
  no_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
            '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23',
            'X', 'Y']
  assert list(utils.get_chromosomes(prepend='')) == no_chr

  # Test with default prepend
  with_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
              'chr23', 'chrX', 'chrY']
  assert list(utils.get_chromosomes()) == with_chr

  # Test including custom contigs
  custom = ['robins_contig', 'space_contig2']
  with_custom = with_chr + custom
  assert list(utils.get_chromosomes(custom=custom)) == with_custom


# def test_get_intervals():
#   # Connect to the fixtures database and get all intervals belonging to
#   # contigs (chromosomes) 'chr1', 'chr2'.
#   db_path = 'sqlite:///tests/fixtures/db.sqlite'
#   db = ElementAdapter(db_path)
#   contig_id = 'chr1'
#   answer = [[(10, 20, 1), (15, 30, 2), (45, 51, 3)]]
#   assert list(utils.get_intervals(db, contig_id)) == answer


def test_group_intervals():
  # Random intervals
  intervals = [(5, 15, 1), (5, 20, 2), (25, 30, 3), (30, 65, 4), (60, 82, 5),
               (83, 88, 6)]

  # Test without extension
  answer = [[(5, 15, 1), (5, 20, 2), (25, 30, 3)],
            [(30, 65, 4)],
            [(60, 82, 5), (83, 88, 6)]]
  result = list(utils.group_intervals(intervals, threshold=30))
  assert result == answer

  # Test with extension +/- 2 bases
  answer = [[(3, 17, 1), (3, 22, 2), (23, 32, 3)],
            [(28, 67, 4)],
            [(58, 84, 5)],
            [(81, 90, 6)]]
  result = list(utils.group_intervals(intervals, threshold=30, extension=2))
  assert result == answer

  # This should also test the edge case when extension make the interval
  # go to <= 0. Throw error!


def test_merge_intervals():
  intervals = [(10, 21), (25, 35)]
  assert utils.merge_intervals(intervals) == (10, 35)


def test_calculate_completeness():
  read_depths = np.array([5., 10., 10., 5., 5., 0.])
  assert utils.calculate_completeness(read_depths, threshold=10) == (2 / 6)

  # Test all 0
  read_depths = np.array([0, 0, 0, 0, 0])
  assert utils.calculate_completeness(read_depths, threshold=1) == 0

  # Test empty array
  read_depths = np.array([])
  assert utils.calculate_completeness(read_depths, threshold=1) == 0


def test_assign_relative_positions():
  # Test regular case
  abs_start = 100
  abs_end = 157
  overall_start = 50
  rel_start = abs_start - overall_start
  rel_end = abs_end - overall_start

  answer = (rel_start, rel_end)
  result = utils.assign_relative_positions(abs_start, abs_end, overall_start)

  assert result == answer

  # Test edge case
  abs_start = 0
  abs_end = 78
  overall_start = 0
  rel_start = abs_start - overall_start
  rel_end = abs_end - overall_start

  answer = (rel_start, rel_end)
  result = utils.assign_relative_positions(abs_start, abs_end, overall_start)

  assert result == answer


def test_calculate_values():
  # Test calculating some coverage and completeness
  read_depths = np.array([5, 10, 20, 20, 10])
  cutoff = 10
  coverage = (5 + 10 + 20 + 20 + 10) / 5
  completeness = 4 / 5

  result = utils.calculate_values(read_depths, threshold=cutoff)
  assert result == (coverage, completeness)
