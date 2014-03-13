# -*- coding: utf-8 -*-

import pytest
from .context import CoverageAdapter


bam = CoverageAdapter('tests/fixtures/align.bam')


def rest_read():
  # These are the main results from 'align.bam'
  depths = [2., 4., 5., 5., 5., 5., 6., 7., 7., 7., 7., 7., 7., 7., 7., 7., 7.,
            7., 7., 7., 7., 7., 7., 8., 8., 7., 7., 7., 7., 7., 7., 7., 6., 4.,
            4., 3., 3., 2., 2.]

  # Read BAM from position [1,10]
  depths = bam.read('chr1', 0, 9)

  # Make assertions: we expect the read depths from 1st to 10th pos to be
  # included.
  answer = depths[0:10]
  assert list(depths) == answer

  # Test also an interval that extends beyond the reads
  depths = bam.read('chr1', 35, 45)
  assert list(depths) == depths[35:] + [0, 0, 0, 0, 0, 0, 0]

  # Test interval completely outside the reads
  depths = bam.read('chr1', 50, 60)
  assert list(depths) == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

  # Test submitting a false chromosome ID
  with pytest.raises(ValueError):
    bam.read('crh1', 10, 20)
