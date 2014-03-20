# -*- coding: utf-8 -*-

"""
chanjo.sql.models
~~~~~~~~~~~~~~~~~~
"""

from datetime import datetime

from sqlalchemy import Table, Column, ForeignKey, String, Integer, DateTime
from sqlalchemy import Text, Float, Boolean
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
  """The :class:`Superset` represents a collection of sets and potentially
  overlapping intervals. It can be related to multiple sets and multiple
  intervals.

  Args:
    superset_id (str): Unique superset Id e.g. HGNC gene symbol
    contig (str): Contig/Chromosome Id
    start (int): 1-based start of the superset (first interval)
    end (int): 1-based end of the superset (last interval, no UTR)
    strand (str): Strand +/-
    secondary_id (str): E.g. Entrez gene Id
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
  """The :class:`Set` represents a set of non-overlapping intervals. It can
  *only* be related to a single superset.
  
  Args:
    set_id (str): Unique set Id (e.g. CCDS transcript Id)
    contig_id (str): Contig/Chromosome Id
    start (int): 1-based start of the set (first interval)
    end (int): 1-based end of the set (last interval, no UTR)
    strand (str): Strand +/-
    superset_id (str): Related superset Id, e.g. HGNC gene symbol
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
    """Returns the combined number of exon bases. Excludes intronic bases.

    Returns:
      int: Total 'intervalic' (exonic) length of the set
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

  It can be related to a multiple :class:`Set` (transcripts). Start and end
  coordinates are 1-based.

  Args:
    interval_id (str): Unique interval Id
    contig_id (str): Contig/Chromosome Id
    start (int): 1-based start of the interval
    end (int): 1-based end of the interval
    strand (str): Strand +/-
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
    """<magic> Returns the number of bases.

    Returns:
      int: Length of interval in number of bases
    """
    # We add +1 because we count positions and both coordinates are 1-based
    return (self.end - self.start) + 1


# +--------------------------------------------------------------------+
# | Sample ORM classes
# +--------------------------------------------------------------------+
class Sample(Base):
  """Stores metadata about each sample. This helps out in consolidating
  important information in one place.

  .. versionadded:: 0.4.0

  Args:
    sample_id (str): Unique sample Id
    group_id (str): Unique group Id
    cutoff (int): Cutoff used for completeness
    extension (bool): Number of bases added to each interval
    coverage_source (str): Path to the BAM file used
  """
  __tablename__ = 'sample'

  id = Column(String(32), primary_key=True)
  group_id = Column(String(32))

  cutoff = Column(Integer)
  extension = Column(Integer)
  coverage_source = Column(Text)
  element_source = Column(Text)
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

  Args:
    parent_id (int): Parent record ID
    sample_id (str): Unique sample identifier
    group_id (str): Group identifier
    coverage (float): Average coverage for the exon
    completeness (float): Ratio of adequately covered bases
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

  Args:
    parent_id (int): Set ID
    sample_id (str): Unique sample identifier
    group_id (str): Group identifier
    coverage (float): Average coverage for the set
    completeness (float): Ratio of adequately covered bases
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

  Args:
    parent_id (int): Superset ID
    sample_id (str): Unique sample identifier
    group_id (str): Group identifier
    coverage (float): Average coverage for the superset
    completeness (float): Ratio of adequately covered bases
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
