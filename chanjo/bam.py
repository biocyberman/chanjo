#!/usr/bin/env python
# coding: utf-8
"""
  bam.module
  ~~~~~~~~~~~~~

  Responsible for returning BEDGraph intervals given start, end coordinates on
  a specified chromosome.

  :copyright: 2013 by Robin Andeer, see AUTHORS for more details
  :license: license_name, see LICENSE for more details
"""

import pysam


class CoverageAdaptor(pysam.Samfile):
  """
  Chanjo adaptor for fetching BEDGraph intervals directly from BAM alignment
  files.

  :param bam_path: [str] Path to the BAM alignment file. This is required at
                   setup.

  Usage:
    from chanjo.bam import CoverageAdaptor
    path = "/path/to/bam/file.bam"
    adaptor = CoverageAdaptor(path)
  """

  def __init__(self, bam_path):
    super(CoverageAdaptor, self).__init__(bam_path, "rb")

  def coverage(self, chrom, intervals, maxDepth=float("inf"), sort=False):
    """
    Scans all intervals and calculates the combined coverage statistics.
    Expects intervals sorted by start position.
    """

    try:
      # Start positions are 0 based in the genomics world
      firstPos = intervals[0].start - 1
      endPos = intervals[-1].end
    except ValueError:
      return (0, 0)

    counter = Counter(intervals)

    # An empty string is evaluated larger than any integer
    iterator = self.pileup(str(chrom), firstPos, endPos, callback=counter)

    baseCount = float(sum((interval.end - interval.start
                           for interval in intervals)))

    return counter.readCount / baseCount, counter.passedCount / baseCount

  def intervals(self, chrom, intervals, maxDepth=float("inf")):
    """
    Public: Generates BEDGraph intervals of equal coverage between start and end
    on the given chromosome.

    :param chrom:    [str] The chromosome of interest
    :param start:    [int] The first position of the interval
    :param end:      [int] The last position of the interval
    :param maxDepth: [int] The highest read depth to consider

    Usage:
      adaptor.intervals("17", 100023, 102051, 50)
      [out] => [<chanjo.bam.Interval instance at 0x10f2ea518>,
                <chanjo.bam.Interval instance at 0x10f2ea4d0>]
    """
    # Start positions are 0 based in the genomics world
    firstPos = intervals[0].start - 1
    endPos = intervals[-1].end

    # Preallocate an array with enough space (worst case scenario)
    self._intervals = [None]*(endPos-firstPos)
    # Pointer for the array
    self.count = 0

    # An empty string is evaluated larger than any integer
    self.iterator = self.pileup(str(chrom), firstPos, endPos)

    for interval in intervals:
      # Move the iterator to the start of the interval
      lastStart, lastDepth = self._move2start(interval.start)

      # Pick up the iterator
      for col in self.iterator:

        # Don't move beyond the given interval
        if col.pos >= interval.end:
          # Move on to the next interval
          break

        # Testing whether coverage has changed or
        # if both the beginning and the current position in the current interval
        # is above the specified max.
        if col.n != lastDepth and (col.n < maxDepth or lastDepth < maxDepth):

          # Stuff the latest interval into the array
          self._persistCoverage(chrom, lastStart, col.pos, lastDepth, maxDepth)

          # Start new interval
          lastDepth = col.n
          lastStart = col.pos

          # Move the save pointer one step forward
          self.count += 1

      # Stuff the last interval into the array
      self._persistCoverage(chrom, lastStart, col.pos, lastDepth, maxDepth)
      # `count` now is the number of intervals that was created
      self.count += 1

    # Return the subset of the list actually containing calculated intervals
    return self._intervals[:self.count]

  def _move2start(self, start):
    # Move the iterator until start of given interval
    for col in self.iterator:
      if col.pos < start:
        continue
      else:
        # Return the position and read depth
        return col.pos, col.n

  def intervals_old(self, chrom, start, end, maxDepth=float("inf")):
    """
    Public: Generates BEDGraph intervals of equal coverage between start and end
    on the given chromosome.

    :param chrom:    [str] The chromosome of interest
    :param start:    [int] The first position of the interval
    :param end:      [int] The last position of the interval
    :param maxDepth: [int] The highest read depth to consider

    Usage:
      adaptor.intervals("17", 100023, 102051, 50)
      [out] => [<chanjo.bam.Interval instance at 0x10f2ea518>,
                <chanjo.bam.Interval instance at 0x10f2ea4d0>]
    """

    # Start positions are 0 based in the genomics world
    start -= 1

    # An empty string is evaluated larger than any integer
    iterator = self.pileup(str(chrom), start, end)

    lastDepth = None
    lastStart = None

    # Move the iterator until start of given interval
    for col in iterator:
      if col.pos < start:
        continue
      else:
        # Initialize using the first position in the interval
        lastDepth = col.n
        lastStart = col.pos

        break

    if not lastDepth:
      # Pilup came up empty for the interval
      return []

    # Preallocate an array and hope we have enough space
    self._intervals = [None]*50*(end - start)
    # Pointer for the array
    self.count = 0

    # Pick up the iterator and go through one position at a time
    for col in iterator:

      # Don't move beyond the given interval
      if col.pos >= end:
        break

      # Testing whether coverage has changed or
      # if both the beginning and the current position in the current interval
      # is above the specified max.
      if col.n != lastDepth and (col.n < maxDepth or lastDepth < maxDepth):

        # Stuff the latest interval into the array
        self._persistCoverage(chrom, lastStart, col.pos, lastDepth, maxDepth)

        # Start new interval
        lastDepth = col.n
        lastStart = col.pos

        # Move the save pointer one step forward
        self.count += 1

    # Stuff the last interval into the array
    self._persistCoverage(chrom, lastStart, col.pos, lastDepth, maxDepth)
    # `count` now is the number of intervals that was created
    self.count += 1

    # Return the subset of the list actually containing calculated intervals
    return self._intervals[:self.count]

  def _persistCoverage(self, chrom, lastStart, currentPos, lastDepth, maxDepth):
    """
    Private.
    """

    if lastDepth > 0:

      # Limit the max reported read depth
      if lastDepth > maxDepth:
        reportedDepth = maxDepth
      else:
        reportedDepth = lastDepth

      self._intervals[self.count] = Interval(lastStart, currentPos,
                                             reportedDepth)
    else:
      print("Positions with 0 reads: {}".format(currentPos))

class Interval(object):
  """docstring for Interval"""
  def __init__(self, start, end, value=None, chrom=None):
    super(Interval, self).__init__()
    self.start = start
    self.end = end
    self.value = value

  def __len__(self):
    return self.end - self.start

  def __str__(self):
    return "({0} - {1}]".format(self.start, self.end)

class Counter(object):
  def __init__(self, intervals, maxDepth=50):
    # Initialize
    self.intervals = intervals
    self.intervalCount = 0
    self.readCount = 0
    self.passedCount = 0
    self.maxDepth = maxDepth
    self.endPos = intervals[-1].end

    # BEDGraph interval helper variables: keeps track of the current bgInterval
    self.lastStart = None
    self.lastDepth = None

  @property
  def currentInterval(self):
    return self.intervals[self.intervalCount]

  def annotate(self, currPos):
    bases = currPos - self.lastStart

    self.readCount += self.lastDepth * bases

    # If sufficient read depth
    if self.lastDepth >= self.maxDepth:
      self.passedCount += bases

  def __call__(self, col):
    if col.pos <= self.endPos:
      # If the iterator has entered the current interval
      if col.pos >= self.currentInterval.start:

        if self.lastStart is None:
          self.lastStart = col.pos
          self.lastDepth = col.n

        # If the iterator has passed the current interval
        if col.pos > self.currentInterval.end:

          self.intervalCount += 1

          # End the bgInterval prematurely
          self.annotate(col.pos)

          # Update the helper variables
          self.lastDepth = col.n
          self.lastStart = col.pos

        elif col.n != self.lastDepth:

          # The bgInterval ended on the previous position
          self.annotate(col.pos)

          # Update the helper variables
          self.lastDepth = col.n
          self.lastStart = col.pos
