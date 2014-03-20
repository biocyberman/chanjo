# -*- coding: utf-8 -*-

"""
chanjo.bam
~~~~~~~~~~~

The default :class:`CoverageAdapter` that ships with Chanjo. Talks directly
to a BAM alignment file to extract read depth data.

Depends on the Pysam_ package and requires Samtools_ to be installed in your
``$PATH`` to work. N.B. by installing *Pysam* through *pip*, *Samtools* will
by installed alongside.

:copyright: (c) 2013 by Robin Andeer
:license: MIT, see LICENSE for more details

.. _Pysam: http://www.cgat.org/~andreas/documentation/pysam/contents.html
.. _Samtools: http://samtools.sourceforge.net/
"""

import numpy as np
from path import path
from pysam import Samfile


class CoverageAdapter(Samfile):
  """
  Adapter for interfacing directly with BAM alignment files. Inherits from
  :class:`Samfile` which requires a second init parameter that tells it to
  expect a *binary* BAM-file rather than the plain text cousin SAM-file format.

  .. code-block:: python

    >>> from chanjo.bam import CoverageAdapter
    >>> bam_path = "/path/to/bam/file.bam"
    >>> adapter = CoverageAdapter(bam_path)

  Args:
    bam_path (str): Path to a BAM alignment file
  """
  def __init__(self, bam_path):
    super(CoverageAdapter, self).__init__()

    # Raise an error if the file doesn't exist
    if not path(bam_path).exists():
      raise OSError(errno.ENOENT, bam_path)

  def read(self, contig_id, start, end):
    """Generates a list of read depths for each position between (start, end).
    The `numpy` array is used to optimize performance when building and
    slicing the list.

    This method depends on `Pysam` >=0.7.5 since the `truncate` option wasn't
    available in previous versions.

    .. code-block:: python

      >>> adapter.read('17', 0, 5)
      array([3., 3., 4., 4., 5., 4.])

    .. note::

      Positions are expected to be 0:0-based. In other words; if start=0,
      end=9 you should expect read depths for base pair positions 1-10 to
      be returned.

    Args:
      contig_id (str): The contig/chromosome Id (str) of interest
      start (int): The first position of the interval (0-based)
      end (int): The last position of the interval (0-based)

    Returns:
      numpy.array: Array of read depths for *each* position in the interval
    """
    # Generate a list of 0 read depth for each position
    positions = np.zeros(end + 1 - start)

    # Start Pileup iterator and walk through each position in the interval
    # `truncate` will make sure it starts and ends on the given positions!
    # +1 to end because pysam otherwise stops one base short by default
    for col in self.pileup(str(contig_id), start, end + 1, truncate=True):

      # Overwrite the read depth in the correct position
      # This will allow simple slicing to get at the positions of interest
      positions[col.pos - start] = col.n

    return positions
