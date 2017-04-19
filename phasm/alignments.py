import enum
from abc import ABCMeta, abstractmethod, abstractproperty
from typing import List, Tuple, NamedTuple


class Strand(enum.IntEnum):
    SAME = 0
    OPPOSITE = 1


class AlignmentType(enum.IntEnum):
    OVERLAP_AB = 0
    OVERLAP_BA = 1
    A_CONTAINED = 2
    B_CONTAINED = 3


class DNASegment(metaclass=ABCMeta):
    """Represents a piece of DNA. Does not necessarily contains the actual
    sequence, its length and identifier is more important.

    This abstract class serves as base classs for more specific types of DNA
    segments."""

    @abstractmethod
    def with_orientation(self, orientation: str) -> 'OrientedDNASegment':
        """Return an `OrientedDNASegment` which includes the orientation of this
        segment."""

    @abstractmethod
    def __len__(self) -> int:
        """Length of this DNA segment"""

    @abstractmethod
    def __str__(self) -> str:
        """Identifier of this DNA segment"""

    def __hash__(self) -> int:
        """Hash for this DNA segment to make sure these objects can be used as
        nodes in a networkx graph."""

        return hash(str(self))


class OrientedDNASegment(DNASegment):
    @abstractproperty
    def orientation(self):
        """Return the orientation of a read. Denoted as a '+' or a '-'."""

    @abstractmethod
    def reverse(self) -> 'OrientedDNASegment':
        """Return the reverse complement of this segment.

        Not the actual reverse complement (because these objects do not include
        the DNA sequence itself), but an identifier denoting the reverse
        complement of this segment."""

    def __hash__(self):
        return hash(str(self))


_Read = NamedTuple(
    'Read',
    [('id', str), ('moviename', str), ('well', str), ('pulse_start', int),
     ('pulse_end', int), ('length', int)]
)


class Read(DNASegment, _Read):
    """Represents a read in the database. Can be used as node in an
    :class:`AssemblyGraph`."""

    def __len__(self) -> int:
        return self.length

    def __str__(self) -> str:
        return self.id

    def with_orientation(self, orientation) -> 'OrientedRead':
        return OrientedRead(self, orientation)


_OrientedRead = NamedTuple('OrientedRead',
                           [('read', Read), ('strand', str)])


class OrientedRead(OrientedDNASegment, _OrientedRead):
    """Represents a read with orientation (from which strand). Can be used as
    node in an :class:`AssemblyGraph`."""

    @property
    def orientation(self) -> str:
        return self.strand

    def reverse(self) -> 'OrientedRead':
        return OrientedRead(
            self.read, "-" if self.orientation == "+" else "+")

    def __getattr__(self, key):
        data = self.read._asdict()

        if key in data:
            return data[key]
        else:
            raise AttributeError("OrientedRead has no attribute '{}'".format(
                key))

    def __len__(self):
        return len(self.read)

    def __str__(self):
        return str(self.read) + self.orientation


_MergedReads = NamedTuple(
    'MergedReads',
    [('id', str), ('length', int), ('strand', str),
     ('reads', List[OrientedRead])]
)


class MergedReads(OrientedDNASegment, _MergedReads):
    """Represents a collection of reads to be merged together. For example,
    this is being used in an :class:`AssemblyGraph` to merge unambiguous
    (non-branching) paths to a single node."""

    @property
    def orientation(self) -> str:
        return self.strand

    def reverse(self) -> 'MergedReads':
        new_strand = "+" if self.strand == "-" else "-"
        return MergedReads(self.id, self.length, new_strand,
                           [r.reverse() for r in reversed(self.reads)])

    def __len__(self):
        return self.length

    def __str__(self):
        return self.id + self.strand


_LocalAlignment = NamedTuple(
    '_LocalAlignment',
    [('a', Read), ('b', Read), ('strand', Strand), ('arange', Tuple[int, int]),
     ('brange', Tuple[int, int]), ('differences', int),
     ('tracepoints', List[Tuple[int, int]])]
)


class LocalAlignment(_LocalAlignment):
    def get_overlap_length(self) -> int:
        return max(self.arange[1] - self.arange[0],
                   self.brange[1] - self.brange[0])

    def get_error_rate(self) -> float:
        return self.differences / self.get_overlap_length()

    def get_overhang(self) -> int:
        return (min(self.arange[0], self.brange[0]) +
                min(len(self.a) - self.arange[1],
                    len(self.b) - self.brange[1]))

    def classify(self) -> AlignmentType:
        if (self.arange[0] <= self.brange[0] and len(self.a)-self.arange[1]
                <= len(self.b)-self.brange[1]):
            return AlignmentType.A_CONTAINED
        elif (self.arange[0] >= self.brange[0] and len(self.a)-self.arange[1]
              >= len(self.b)-self.brange[1]):
            return AlignmentType.B_CONTAINED
        elif self.arange[0] >= self.brange[0]:
            return AlignmentType.OVERLAP_AB
        else:
            return AlignmentType.OVERLAP_BA

    def get_oriented_reads(self) -> Tuple[OrientedRead, OrientedRead]:
        return self.a.with_orientation('+'), self.b.with_orientation(
            '+' if self.strand == Strand.SAME else '-')

    def __len__(self):
        return self.get_overlap_length()
