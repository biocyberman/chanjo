# -*- coding: utf-8 -*-

"""
chanjo.sql.core
~~~~~~~~~~~~~~~~~~
The Central Database Class
"""

from sqlalchemy import create_engine, func
from sqlalchemy.orm import sessionmaker

from .models import (Base, Interval, Set, Superset, Interval_Set,
                     IntervalData, SetData, SupersetData, Sample)


class ElementAdapter(object):
  """:class:`ElementalDB` collects functionality needed to setup and interact
  with a SQLAlchemy session. This is likely the only element user will directly
  interact with.

  .. code-block:: python

    >>> from elemental.core import ElementalDB
    >>> db = ElementalDB('data/elements.sqlite')

  .. note::

    For testing pourposes; use ":memory:" as the `path` argument to set up
    in-memory version of the database.

  :param str path: Path to the database to connect to
  :param bool debug: Whether to print logging information (optional)
  :param str dialect: The type of database to connect to (optional)
-------------------------------------------------------------------------------
  SQLAlchemy-based :class:`ElementAdapter` for Chanjo.

  Inherits the basic SQL table structure from `ElementalDB`.

  .. note::

    For testing pourposes; use ":memory:" as the ``path`` parameter to set up
    an in-memory instance of the database.

  :param str path: Path to the database to connect to
  :param str dialect: Type of database to use: 'sqlite' or 'mysql' (optional)
  :param bool debug: Whether to print logging information (optional)
  """
  def __init__(self, path, dialect='sqlite', debug=False):
    super(ElementAdapter, self).__init__()
    # Save store pointer
    self.path = path

    # Connect tose a database
    if dialect == 'sqlite':
      self.engine = create_engine('sqlite:///' + path, echo=debug)
    else:
      # Build path containing username:password@server/database
      auth_path = '{type}://{path}'.format(type=dialect, path=path)
      self.engine = create_engine(auth_path, pool_recycle=3600, echo=debug)

    # Make sure the same engine is propagated to the Base classes
    Base.metadata.bind = self.engine

    # Start a session
    self.session = sessionmaker(bind=self.engine)()
    # Shortcut to query
    self.query = self.session.query

    # ORM class shortcuts
    self.classes = {
      'superset': Superset,
      'set': Set,
      'interval': Interval,
      'interval_set': Interval_Set,
      'superset_data': SupersetData,
      'set_data': SetData,
      'interval_data': IntervalData,
      'sample': Sample
    }

  def setup(self):
    """<public> Sets up a new database with the default tables and columns.

    :returns: ``self``
    """
    # Create the tables
    Base.metadata.create_all(self.engine)

    return self

  def tare_down(self):
    """<public> Tares down a database (tables and columns).

    :returns: ``self``
    """
    # Create the tables
    Base.metadata.drop_all(self.engine)

    return self

  def get(self, typ, type_id):
    """<public> Fetches a specific element or ORM class. Calls itself
    recursively when asked to fetch an element.

    .. code-block:: python

      >>> db = ElementalDB('path/to/element.sqlite')

      # Get a specific gene from the database
      >>> gene = db.get('gene', 'GIT1')

    :param str typ: Element key or 'class'
    :param str type_id: Element ID or ORM class ID
    :returns: An element or ORM class
    """
    if typ == 'class':
      return self.classes[type_id]

    # Get an ORM class
    klass = self.get('class', typ)

    # Return the requested element object (or ``None``) if not found
    return self.session.query(klass).get(type_id)

  def find(self, klass_id, query=None, attrs=None):
    """<public> If the 'query' parameter is a string `find` will fetch one
    element; just like `get`. If query is a list it will match element IDs to
    items in that list and return a list of elements. If 'query' is ``None``
    all elements of that class will be returned.

    .. versionchanged: 0.2.0

    :param str klass_id: The type of element to find
    :param str/list query: (optional) Element ID(s)
    :param list attrs: (optional) List of columns to fetch
    :returns: Element(s) from the database
    :rtype: object/list
    """
    # Get an ORM class
    klass = self.get('class', klass_id)

    if attrs is not None:
      params = [getattr(klass, attr) for attr in attrs]
    else:
      params = (klass,)

    if query is None:
      # Return all `klass_id` elements in the database
      return self.session.query(*params).all()

    elif isinstance(query, list):
      # Return all `klass_id` elements in the database
      return self.session.query(*params).filter(klass.id.in_(query)).all()

    elif isinstance(query, str):
      # Call 'get' to return the single element
      return self.get(klass_id, query)

    else:
      # For now
      return []

  def add(self, elements):
    """<public> Add one or multiple new elements to the database and commit the
    changes. Chainable.

    :param elements: New ORM object instance or list of such
    :type elements: object or list
    :returns: ``self`` for chainability
    """
    if isinstance(elements, Base):
      # Add the record to the session object
      self.session.add(elements)

    elif isinstance(elements, list):
      # Add all records to the session object
      self.session.add_all(elements)

    return self

  def create(self, class_id, *args, **kwargs):
    """<public> Creates a new instance of an ORM element object filled in with
    the given `attributes`.

    If attributes is a tuple they must be in the correct order. Supplying a
    `dict` doesn't require the attributes to be in any particular order.

    :param str class_id: Choice between "superset", "set", "interval"
    :param \*args: List the element attributes in the *correct order*
    :param \**kwargs: Element attributes in whatever order you like
    :returns: The new ORM instance object
    """
    if args:
      # Unpack tuple
      return self.get('class', class_id)(*args)
    elif kwargs:
      # Unpack dictionary
      return self.get('class', class_id)(**kwargs)
    else:
      raise TypeError('Submit attributes as arguments or keyword arguments')

  def commit(self):
    """<public> Manually persist changes made to various elements. Chainable.

    :returns: ``self`` for chainability
    """
    # Commit/persist dirty changes to the database
    self.session.commit()

    return self

  def set_stats(self, sample_id):
    """<public> Calculates set level metrics to annotate transcripts.
    Requires all related exons to already be properly annotated.

    What's happening is that we are summing read depths and adequately covered
    bases for each exon in a transcript and then dividing those numbers by the
    total exon length of the transcript.

    .. note::

      Transcript annotation needs to be carried out before annotating genes!

    :param str sample_id: Sample ID to match with coverage annotations
    :returns: List of tuples: ``(<tx_id>, <coverage>, <completeness>)``
    """
    # Length of interval (in number of bases, hence +1)
    interval_length = Interval.end - Interval.start + 1
    # Length of interval times mean coverage
    cum_interval_coverage = interval_length * IntervalData.coverage
    # Summed cumulative coverage for all intervals (of the set)
    cum_set_coverage = func.sum(cum_interval_coverage)
    # Summed interval length of all intervals (of the set)
    total_set_length = func.sum(interval_length)
    # Cumulative set coverage divided by total set length
    mean_set_coverage = cum_set_coverage / total_set_length

    # Length of interval times mean completeness
    cum_interval_completeness = interval_length * IntervalData.completeness
    # Summed cumulative completeness for all intervals (of the set)
    cum_set_completeness = func.sum(cum_interval_completeness)
    # Cumulative set coverage divided by total set length
    mean_set_completeness = cum_set_completeness / total_set_length

    # Values to fetch
    interval_set_columns = Interval_Set.columns.values()
    interval_id = interval_set_columns[0]
    set_id = interval_set_columns[1]

    return self.query(
      set_id,
      mean_set_coverage,
      mean_set_completeness
    ).join(Interval, interval_id == Interval.id)\
     .join(IntervalData, Interval.id == IntervalData.parent_id)\
     .filter(IntervalData.sample_id == sample_id).group_by(set_id)

  def superset_stats(self, sample_id):
    """<public> Calculates superset level metrics to annotate genes.
    Requires all related sets to already be properly annotated.

    What's happening is that we are simply taking the average of the metrics
    on the transcript level and applying that as gene metrics. This gives a
    decent, albeit not perfect, represenation of gene level metrics.

    .. note::

      Annotation of transcripts needs to be acomplished before annotating
      genes!

    :param str sample_id: Sample ID to match with coverage annotations
    :returns: List of tuples: ``(<gene_id>, <coverage>, <completeness>)``
    """
    return self.query(
      Set.superset_id,
      func.avg(SetData.coverage),
      func.avg(SetData.completeness)
    ).join(SetData, Set.id == SetData.parent_id)\
     .filter(SetData.sample_id == sample_id)\
     .group_by(Set.superset_id)
