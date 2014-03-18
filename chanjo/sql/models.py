# -*- coding: utf-8 -*-

"""
chanjo.sql.models
~~~~~~~~~~~~~~~~~~
"""

from datetime import datetime

from sqlalchemy import (Table, Column, ForeignKey, String, Integer, DateTime,
                        Float, Boolean)
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base

# Base for declaring a mapping
Base = declarative_base()


# +--------------------------------------------------------------------+
# | Association tables
# | ~~~~~~~~~~~~~~~~~~~
# | Provides the many-to-many relationships between:
# | - Interval<->Set
# +--------------------------------------------------------------------+
Interval_Set = Table('interval__set', Base.metadata,
  Column('interval_id', String(32), ForeignKey('interval.id')),
  Column('set_id', String(32), ForeignKey('set.id'))
)


# +--------------------------------------------------------------------+
# | Superset ('gene') ORM
# +--------------------------------------------------------------------+
class Superset(Base):
  """
  The :class:`Superset` represents a collection of sets and potentially
  overlapping intervals (gene). It can be related to multiple sets and
  multiple intervals.

  :param str superset_id: E.g. HGNC gene symbol
  :param str contig: Contig/Chromosome ID
  :param int start: 0-based start of the gene (first interval)
  :param int end: 0-based end of the gene (last interval, no UTR)
  :param str strand: Strand +/-
  :param str secondary_id: E.g. Entrez gene ID
  """
  __tablename__ = 'superset'

  id = Column(String(32), primary_key=True)
  secondary_id = Column(Integer)
  contig_id = Column(String(5))
  start = Column(Integer)
  end = Column(Integer)
  strand = Column(String(1))

  def __init__(self, superset_id, contig_id, start, end, strand, secondary_id):
    super(Superset, self).__init__()

    self.id = superset_id
    self.secondary_id = secondary_id
    self.contig_id = contig_id
    self.start = start
    self.end = end
    self.strand = strand


# +--------------------------------------------------------------------+
# | Set ('transcript') ORM
# +--------------------------------------------------------------------+
class Set(Base):
  """The :class:`Set` represents a set of non-overlapping intervals
  (transcript). It can *only* be related to a single superset.

  :param str set_id: Unique set ID (e.g. CCDS transcript ID)
  :param str contig_id: Contig/Chromosome ID
  :param int start: 0-based start of the transcript (first interval)
  :param int end: 0-based end of the transcript (last interval, no UTR)
  :param str strand: Strand +/-
  :param str superset_id: E.g. HGNC gene symbol
  """
  __tablename__ = 'set'

  id = Column(String(32), primary_key=True)
  contig_id = Column(String(5))
  start = Column(Integer)
  end = Column(Integer)
  strand = Column(String(1))

  superset_id = Column(String(32), ForeignKey('superset.id'))
  superset = relationship(Superset, backref=backref('sets', order_by=start))

  def __init__(self, set_id, contig_id, start, end, strand, superset_id=None):
    super(Set, self).__init__()

    self.id = set_id
    self.contig_id = contig_id
    self.start = start
    self.end = end
    self.strand = strand
    self.superset_id = superset_id

  def __len__(self):
    """<magic> Returns the combined number of exon bases. Excludes intronic
    bases.
    """
    base_count = 0
    for interval in self.intervals:
      base_count += len(interval)

    return base_count


# +--------------------------------------------------------------------+
# | Interval ('exon') ORM
# +--------------------------------------------------------------------+
class Interval(Base):
  """The :class:`Interval` represents a continous genetic interval on a given
  contig (e.g. exon).

  It can be related to a multiple :class:`Set` (transcripts) and multiple
  :class:`Superset` (genes). Start and end coordinates are 0-based.

  :param str interval_id: Unique interval ID
  :param str contig_id: Contig/Chromosome ID
  :param int start: 0-based start of the interval
  :param int end: 0-based end of the interval
  :param str strand: Strand +/-
  """
  __tablename__ = 'interval'

  id = Column(String(32), primary_key=True)
  contig_id = Column(String(5))
  start = Column(Integer)
  end = Column(Integer)
  strand = Column(String(1))
  first = Column(Boolean)
  last = Column(Boolean)

  # This also defines the ``backref`` to give transcripts an exons property
  sets = relationship(Set, secondary=Interval_Set,
                      backref=backref('intervals', order_by=start))

  def __init__(self, interval_id, contig_id, start, end, strand):
    super(Interval, self).__init__()

    self.id = interval_id
    self.contig_id = contig_id
    self.start = start
    self.end = end
    self.strand = strand

  def __len__(self):
    """
    <magic> Returns the number of bases.
    """
    # We add +1 because we count positions and both coordinates are 0-based
    return (self.end - self.start) + 1


# +--------------------------------------------------------------------+
# | Sample ORM classes
# +--------------------------------------------------------------------+
class Sample(Base):
  """Stores metadata about each sample. This helps out in consolidating
  important information in one place.

  .. versionadded:: 0.4.0

  :param str sample_id: Unique sample ID
  :param str group_id: Unique group ID
  :param int cutoff: Cutoff used for completeness
  :param bool extension: Number of bases added to each interval
  :param str coverage_source: Path to the BAM file used
  """
  __tablename__ = 'sample'

  id = Column(String(32), primary_key=True)
  group_id = Column(String(32))

  cutoff = Column(Integer)
  extension = Column(Integer)
  coverage_source = Column(String(32))
  element_source = Column(String(32))
  created_at = Column(DateTime, default=datetime.now)
  updated_at = Column(DateTime, default=datetime.now, onupdate=datetime.now)

  def __init__(self, sample_id, group_id, cutoff, coverage_source,
               element_source, extension):
    super(Sample, self).__init__()

    self.id = sample_id
    self.group_id = group_id
    self.cutoff = cutoff
    self.coverage_source = coverage_source
    self.element_source = element_source
    self.extension = extension


# +--------------------------------------------------------------------+
# | Interval Data ORM
# +--------------------------------------------------------------------+
class IntervalData(Base):
  """Stores coverage metrics for a single interval and a given sample. It has a
  many-to-one relationship with it's parent interval object through it's
  ``parent_id`` attribute.

  :param int parent_id: Parent record ID
  :param str sample_id: Unique sample identifier
  :param str group_id: Group identifier
  :param float coverage: Average coverage for the exon
  :param float completeness: Ratio of adequately covered bases
  """
  __tablename__ = 'interval_data'

  id = Column(Integer, primary_key=True, autoincrement=True)
  coverage = Column(Float)
  completeness = Column(Float)

  # These columns map coverage/completeness to sample+group
  sample_id = Column(String(32), ForeignKey('sample.id'))
  sample = relationship(Sample, backref=backref('intervals'))
  group_id = Column(String(32))

  # Genetic relationship
  parent_id = Column(String(32), ForeignKey('interval.id'))
  parent = relationship(Interval, backref=backref('data'))

  def __init__(self, parent_id, sample_id=None, group_id=None,
               coverage=None, completeness=None):
    super(IntervalData, self).__init__()

    self.parent_id = parent_id
    self.sample_id = sample_id
    self.group_id = group_id
    self.coverage = coverage
    self.completeness = completeness


# +--------------------------------------------------------------------+
# | Set Data ORM
# +--------------------------------------------------------------------+
class SetData(Base):
  """Stores coverage metrics for a single set and a given sample. It has a
  many-to-one relationship with it's parent set object through it's
  ``parent_id`` attribute.

  :param int parent_id: Set ID
  :param str sample_id: Unique sample identifier
  :param str group_id: Group identifier
  :param float coverage: Average coverage for the set
  :param float completeness: Ratio of adequately covered bases
  """
  __tablename__ = 'set_data'

  id = Column(Integer, primary_key=True, autoincrement=True)
  coverage = Column(Float)
  completeness = Column(Float)

  # These columns map coverage/completeness to an individual+group
  sample_id = Column(String(32), ForeignKey('sample.id'))
  sample = relationship(Sample, backref=backref('sets'))
  group_id = Column(String(32))

  # Genetic relationship
  parent_id = Column(String(32), ForeignKey('set.id'))
  parent = relationship(Set, backref=backref('data'))

  def __init__(self, parent_id, sample_id=None, group_id=None,
               coverage=None, completeness=None):
    super(SetData, self).__init__()

    self.parent_id = parent_id
    self.sample_id = sample_id
    self.group_id = group_id
    self.coverage = coverage
    self.completeness = completeness


# +--------------------------------------------------------------------+
# | Superset Data ORM classes
# +--------------------------------------------------------------------+
class SupersetData(Base):
  """Stores coverage metrics for a single superset and a given sample. It has a
  many-to-one relationship with it's parent superset object through it's
  ``parent_id`` attribute.

  :param int parent_id: Superset ID
  :param str sample_id: Unique sample identifier
  :param str group_id: Group identifier
  :param float coverage: Average coverage for the superset
  :param float completeness: Ratio of adequately covered bases
  """
  __tablename__ = 'superset_data'

  id = Column(Integer, primary_key=True, autoincrement=True)
  coverage = Column(Float)
  completeness = Column(Float)

  # These columns map coverage/completeness to an individual+group
  sample_id = Column(String(32), ForeignKey('sample.id'))
  sample = relationship(Sample, backref=backref('supersets'))
  group_id = Column(String(32))

  # Genetic relationship
  parent_id = Column(String(32), ForeignKey('superset.id'))
  parent = relationship(Superset, backref=backref('data'))

  def __init__(self, parent_id, sample_id=None, group_id=None,
               coverage=None, completeness=None):
    super(SupersetData, self).__init__()

    self.parent_id = parent_id
    self.sample_id = sample_id
    self.group_id = group_id
    self.coverage = coverage
    self.completeness = completeness
