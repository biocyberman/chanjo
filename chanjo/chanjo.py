#!/usr/bin/env python
# coding: utf-8
"""
  chanjo.module
  ~~~~~~~~~~~~~

  The clue and control hub of the `chanjo` package.
  Borrows many structural ideas from Ember.js in the way you add in adaptors to
  add in modular functionality that can be switched out to support multiple
  "backends".

  :copyright: 2013 by Robin Andeer, see AUTHORS for more details
  :license: MIT, see LICENSE for more details
"""

from __future__ import print_function
import collections
from interval import Interval, IntervalSet


class Analyzer(object):
  """docstring for Analyzer"""
  def __init__(self, coverageAdaptor=None, elementAdaptor=None):
    super(Analyzer, self).__init__()

    # Set up the adaptors
    if coverageAdaptor and elementAdaptor:
      self.setAdaptors(coverageAdaptor, elementAdaptor)

  def setAdaptors(self, coverageAdaptor, elementAdaptor):
    """
    Public: Plugs in the required adaptors and sets up a few shortcuts.
    ----------
    :param coverageAdaptor: [object] A class instance of a Coverage Adaptor
    :param elementAdaptor: [object] A class instance of a Element Adaptor

    Usage:
      from chanjo.bigBed import CoverageAdaptor
      from chanjo.sqlite import ElementAdaptor

      analyzer = Analyzer()
      bam_path = "/path/to/file.bam"
      cov_path = "/path/to/database.db"
      analyzer.setAdaptors(CoverageAdaptor(bam_path), ElementAdaptor(cov_path))
    """
    # Customizable adaptors
    self.coverageAdaptor = coverageAdaptor
    self.elementAdaptor = elementAdaptor

    # Shortcut to getting elements by ID
    self.get = self.elementAdaptor.get

    # Shortcut to getting coverage for intervals
    self.intervals = self.coverageAdaptor.intervals

  def annotate(self, elem_class, elem_id, cutoff=50):
    """
    Public: Calculate and persist coverage for a single element.
    ----------
    :param elem_class: [str] The class of element, e.g. "gene"
    :param elem_id:    [str] The unique element id
    :param cutoff:     [int] The read depth level to use for coverage
                       completeness (Default: 50)

    Usage:
      analyzer.annotate("gene", "EGFR", 50)
      [out] => <chanjo.sqlite2.Gene at 0x1041c4c10>
    """
    # For some reason, unicode string doesn't work
    element = self.get(elem_class, str(elem_id))
    if elem_class == "gene" or elem_class == "transcript":
      cov, comp = self.coverage(element.chrom, element.simpleIntervals(),
                                cutoff)
    else:
      # Exon only a single interval
      cov, comp = self.coverage(element.chrom, (element,), cutoff)

    # Update the element with the calculated coverage information
    element.coverage = cov
    element.completeness = comp
    element.cutoff = cutoff

    # Persist the changes
    element.save()

    return element

  def coverage(self, chrom, intervals, cutoff=50, bgIntervals=None):
    """
    Doesn't handle overlapping intervals
    """
    # Initialize
    totBaseCount = 0
    readCount = 0
    passedCount = 0

    # Make sure at least one interval was submitted
    try:
      firstPos = intervals[0].start
      lastPos = intervals[-1].end
    except IndexError:
      # No intervals were submitted
      return (0, 0)

    if bgIntervals is None:
      # Get BEDGraph intervals covering all input intervals
      # Minimizes the number of times we have to fetch BEDGraph intervals
      bgIntervals = self.intervals(chrom, intervals)

    bgCount = 1

    for interval in intervals:
      # We need the number of bases and we are using 1-based positions
      # This should be defined universally as either 0,1-based or 1,1-based.
      totBaseCount += interval.end - (interval.start-1)

      # Two BEDGraph intervals could overlap
      bgCount -= 1

      # Continue until we moved passed the current BEDGraph interval
      while bgCount < len(bgIntervals) and bgIntervals[bgCount].start <= interval.end:

        # Test overlap between BEDGraph interval and current input interval
        if interval.end < bgIntervals[bgCount].end:
          bgBases = interval.end - bgIntervals[bgCount].start
        else:
          bgBases = bgIntervals[bgCount].end - bgIntervals[bgCount].start

        # Add the number of overlapping reads
        readCount += (bgBases * bgIntervals[bgCount].value)

        # Add the position if it passes `cutoff`
        if bgIntervals[bgCount].value >= cutoff:
          passedCount += bgBases

        bgCount += 1

    return (readCount / float(totBaseCount),
            passedCount / float(totBaseCount))

  def levels(self, element, intervals=None, save=True):
    """
    Save coverage as discrete levels to be vizualized across genes.
    """
    # It's possible to "reuse" BEDGraph covereage intrevals.
    if not intervals:
      # Assume gene
      intervals = self.getCoverage(element.chrom,
                                   element.intervals.lower_bound(),
                                   element.intervals.upper_bound())

    # Four coverage levels are used, `IntervalSet` will automatically merge
    # the calculated level intervals.
    levels = [IntervalSet(), IntervalSet(), IntervalSet(), IntervalSet()]
    for interval in intervals:
      readDepth = interval.value

      if readDepth > 50:
        levels[3].add(Interval(interval.start, interval.end))
      elif readDepth > 10:
        levels[2].add(Interval(interval.start, interval.end))
      elif readDepth > 0:
        levels[1].add(Interval(interval.start, interval.end))
      else:
        levels[0].add(Interval(interval.start, interval.end))

    for count, level in enumerate(levels):
      # Stringify the merged intervals to be able to store them in SQL database
      intervals = ",".join(["{0}-{1}".format(ival.lower_bound, ival.upper_bound)
                            for ival in level.intervals])

      # Place the stringified intervals under the correct attribute
      if count == 0:
        element.intervals_0x = intervals
      elif count == 1:
        element.intervals_1_10x = intervals
      elif count == 2:
        element.intervals_10_50x = intervals
      else:
        element.intervals_50x = intervals

    if save:
      # Persist changes
      element.save()

    # Enable chaining
    return element
