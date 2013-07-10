import autumn
from autumn.util import AutoConn
from autumn.model import Model
from autumn.db.relations import ForeignKey, OneToMany
from autumn.db.query import Query
from interval import Interval, IntervalSet
import collections


class ElementAdaptor(object):
  """docstring for ElementAdaptor"""
  def __init__(self, db_path):
    super(ElementAdaptor, self).__init__()

    self.classes = None

    # Set up connection to database
    self.connect(db_path)

  def connect(self, path, new=False):
    # get a database connection object
    # Doesn't have to exist in the first place!
    self.db = AutoConn(path)

    # Unless setting up a new database
    if not new:
      self._defaultORM()

  def get(self, elem_class, elem_ids):
    """
    Get one or multiple elements by ID.
    """
    singleElement = False
    # Test if just a single ID was submitted
    if isinstance(elem_ids, str):
      elem_ids = (elem_ids, )
      singleElement = True
    else:
      # Otherwise make sure we have an iterable
      if not isinstance(elem_ids, collections.Iterable):
        print("Input must be ID string or list of IDs.")
        return None

    elements = range(0,len(elem_ids))
    for count, elem_id in enumerate(elem_ids):

      elements[count] = self.classes[elem_class].get(elem_id)

      # if elements[count]:
      #   # Custom many-to-many relationships
      #   if elem_class == "transcript":
      #     elements[count].exons = self._convertTxEx(transcript_id=elem_id)
      #   elif elem_class == "exon":
      #     elements[count].transcripts = self._convertTxEx(exon_id=elem_id)

    # Return a single element if that was requested
    if singleElement:
      return elements[0]
    else:
      # If a list of element ID was submitted, return a list of elements
      return elements

  def set(self, elem_class, elem_tuple):
    self.classes[elem_class](*elem_tuple).save()

  def setup(self, gs_cols="", tx_cols="", ex_cols=""):
    # code to create the database tsable
    # This is where custom annotations whould be added!
    if len(gs_cols) == 0:
      comma = ""
    else:
      comma = ","

    _genes_sql = """
    DROP TABLE IF EXISTS gene;
    CREATE TABLE gene (
      id TEXT PRIMARY KEY,
      chrom TEXT,
      strand TEXT,
      coverage REAL,
      completeness REAL,
      cutoff INT{0}{1}
    );""".format(comma, gs_cols)

    if len(tx_cols) == 0:
      comma = ""
    else:
      comma = ","

    _transcripts_sql = """
    DROP TABLE IF EXISTS transcripts;
    CREATE TABLE transcripts (
      id TEXT PRIMARY KEY,
      chrom TEXT,
      strand TEXT,
      gene_id TEXT,
      coverage REAL,
      completeness REAL,
      cutoff INT,
      FOREIGN KEY (gene_id) REFERENCES gene(id){0}{1}
    );""".format(comma, tx_cols)
    
    if len(ex_cols) == 0:
      comma = ""
    else:
      comma = ","

    _exons_sql = """
    DROP TABLE IF EXISTS exons;
    CREATE TABLE exons (
      id TEXT PRIMARY KEY,
      chrom TEXT,
      strand TEXT,
      gene_id TEXT,
      start INT,
      end INT,
      coverage REAL,
      completeness REAL,
      cutoff INT,
      FOREIGN KEY (gene_id) REFERENCES gene(id){0}{1}
    );""".format(comma, ex_cols)

    _transcripts_exons_sql = """
    DROP TABLE IF EXISTS transcripts_exons;
    CREATE TABLE transcripts_exons (
      transcript_id TEXT,
      exon_id TEXT,
      PRIMARY KEY (transcript_id, exon_id),
      FOREIGN KEY (transcript_id) REFERENCES transcripts(id),
      FOREIGN KEY (exon_id) REFERENCES exons(id)
    );"""

    # create the tables, dropping any previous tables of the same name
    autumn.util.create_table(self.db, _genes_sql)
    autumn.util.create_table(self.db, _transcripts_sql)
    autumn.util.create_table(self.db, _exons_sql)
    autumn.util.create_table(self.db, _transcripts_exons_sql)

    self._defaultORM()

  def _defaultORM(self):
    # create ORM class; Autumn introspects the database to find out columns
    class Gene(Model):
      db = self.db
      transcripts = OneToMany("Transcript")
      _intervals = None
      _exons = None

      @property
      def exons(self):
        if self._exons is None:
          self._exons = Query(model=Exon).filter(gene_id=self.id).order_by("start")

        return self._exons

      @property
      def exonLength(self):
        """
        Returns the combined length of all related exons.
        N.B. Modifying the __len__ breaks something...
        """
        return sum([interval.upper_bound - interval.lower_bound
                    for interval in self.intervals])

      @property
      def intervals(self):
        """
        Returns all the non-overlapping exonic intervals.
        """
        if not self._intervals:
          self._intervals = IntervalSet([Interval(exon.start, exon.end)
                                        for exon in self.exons])

        return self._intervals

      def simpleIntervals(self):
        return [SimpleInterval(i.lower_bound, i.upper_bound)
                for i in self.intervals]

    class Transcript(Model):
      db = self.db
      gene = ForeignKey(Gene)
      _exons = None

      @property
      def exons(self):
        if not self._exons:
          query = Query(model=Transcript_Exon)
          self._exons = [combo.exon for combo in
                         query.filter(transcript_id=self.id).order_by("start")]

        return self._exons        

      @property
      def exonLength(self):
        """
        Returns the combined length of all related exons.
        """
        return sum([len(exon) for exon in self.exons])

      def simpleIntervals(self):
        """
        Exons aleady have start and end attributes.
        """
        return self.exons

      class Meta:
        table = "transcripts"

    class Exon(Model):
      db = self.db
      gene = ForeignKey(Gene)
      _transcripts = None

      def __len__(self):
        return self.end - self.start

      @property
      def transcripts(self):
        if not self._transcripts:
          query = Query(model=Transcript_Exon)
          self._transcripts = [combo.transcript for combo in
                               query.filter(exon_id=self.id)]

        return self._transcripts

      class Meta:
        table = "exons"

    class Transcript_Exon(Model):
      db = self.db
      transcript = ForeignKey(Transcript, field="transcript_id")
      exon = ForeignKey(Exon, field="exon_id")

      class Meta:
        table = "transcripts_exons"

    # Shortcuts
    self.classes = {
      "gene": Gene,
      "transcript": Transcript,
      "exon": Exon,
      "transcript_exon": Transcript_Exon
    }  

# Very simple interval struct
class SimpleInterval:
  def __init__(self, start, end):
    self.start = start
    self.end = end