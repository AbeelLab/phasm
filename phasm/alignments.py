import enum
from abc import ABCMeta, abstractmethod
from typing import List, Tuple, Union

import numpy

Tracepoints = List[int]
CIGAR = str
Alignment = Union[Tracepoints, CIGAR]
Range = Tuple[int, int]


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

    @abstractmethod
    def __eq__(self, other: 'DNASegment') -> bool:
        pass

    def __hash__(self) -> int:
        """Hash for this DNA segment to make sure these objects can be used as
        nodes in a networkx graph."""

        return hash(str(self))


class OrientedDNASegment(metaclass=ABCMeta):
    @property
    @abstractmethod
    def orientation(self):
        """Return the orientation of a read. Denoted as a '+' or a '-'."""

    @abstractmethod
    def reverse(self) -> 'OrientedDNASegment':
        """Return the reverse complement of this segment.

        Not the actual reverse complement (because these objects do not include
        the DNA sequence itself), but an identifier denoting the reverse
        complement of this segment."""

    @abstractmethod
    def __len__(self) -> int:
        """Length of this DNA segment"""

    @abstractmethod
    def __str__(self) -> str:
        """Identifier of this DNA segment"""

    @abstractmethod
    def __eq__(self, other: 'OrientedDNASegment') -> bool:
        pass

    def __hash__(self) -> int:
        """Hash for this DNA segment to make sure these objects can be used as
        nodes in a networkx graph."""

        return hash(str(self))


class Read(DNASegment):
    """Represents a read in the database. Can be used as node in an
    :class:`AssemblyGraph`."""

    def __init__(self, id: str, length: int, sequence: numpy.array=None):
        self.id = id
        self.length = length
        self.sequence = sequence

    def __len__(self) -> int:
        return self.length

    def __str__(self) -> str:
        return self.id

    def __repr__(self) -> str:
        return "Read[id={0.id}, len={0.length}]".format(self)

    def __eq__(self, other) -> bool:
        if not isinstance(other, self.__class__):
            return False

        return self.id == other.id and self.length == other.length

    def __hash__(self) -> int:
        """Hash for this DNA segment to make sure these objects can be used as
        nodes in a networkx graph."""

        return hash(str(self))

    def with_orientation(self, orientation) -> 'OrientedRead':
        return OrientedRead(self, orientation)


class OrientedRead(OrientedDNASegment):
    """Represents a read with orientation (from which strand). Can be used as
    node in an :class:`AssemblyGraph`."""

    def __init__(self, read: Read, strand: str):
        self.read = read
        self.strand = strand

    @property
    def orientation(self) -> str:
        return self.strand

    def reverse(self) -> 'OrientedRead':
        return OrientedRead(
            self.read, "-" if self.orientation == "+" else "+")

    def __getattr__(self, key):
        try:
            return getattr(self.read, key)
        except AttributeError:
            raise AttributeError("OrientedRead has no attribute '{}'".format(
                key))

    def __len__(self) -> int:
        return len(self.read)

    def __str__(self) -> str:
        return str(self.read) + self.orientation

    def __repr__(self) -> str:
        return "OrientedRead[id={0.id}{0.strand}, len={0.read.length}]".format(
            self)

    def __eq__(self, other) -> bool:
        if not isinstance(other, self.__class__):
            return False

        return self.read == other.read and self.strand == other.strand

    def __hash__(self) -> int:
        """Hash for this DNA segment to make sure these objects can be used as
        nodes in a networkx graph."""

        return hash(str(self))


class MergedReads(OrientedDNASegment):
    """Represents a collection of reads merged together. For example,
    this is being used in an :class:`AssemblyGraph` to merge unambiguous
    (non-branching) paths to a single node."""

    def __init__(self, id: str, length: int, strand: str,
                 reads: List[OrientedRead], prefix_lengths: List[int]):
        self.id = id
        self.length = length
        self.strand = strand
        self.reads = reads

        # Contains a list of integers determining the "unmatched prefix"
        # lengths of each merged read, so we can easily construct the actual
        # DNA sequence of this merged path
        self.prefix_lengths = prefix_lengths

    @property
    def orientation(self) -> str:
        return self.strand

    def reverse(self) -> 'MergedReads':
        new_strand = "+" if self.strand == "-" else "-"

        return MergedReads(self.id, self.length, new_strand,
                           [r.reverse() for r in reversed(self.reads)],
                           list(self.prefix_lengths))

    def __hash__(self) -> int:
        return hash(str(self))

    def __eq__(self, other) -> bool:
        if not isinstance(other, self.__class__):
            return False

        return (self.id == other.id and self.strand == other.strand and
                self.length == other.length)

    def __len__(self) -> int:
        return self.length

    def __str__(self) -> str:
        return self.id + self.strand

    def __repr__(self) -> str:
        return "MergedReads[id={}, r0={}...r{}={}, length={}]".format(
            self.id, self.reads[0], len(self.reads)-1, self.reads[-1],
            self.length
        )


class LocalAlignment:
    def __init__(self, a: Read, b: Read, strand: Strand, arange: Range,
                 brange: Range, alignment: Alignment=None):
        self.a = a
        self.b = b
        self.strand = strand
        self.arange = arange
        self.brange = brange
        self.alignment = alignment

    def get_overlap_length(self) -> int:
        return max(self.arange[1] - self.arange[0],
                   self.brange[1] - self.brange[0])

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

    @property
    def a_id(self) -> str:
        return self.a.id

    @property
    def b_id(self) -> str:
        return self.b.id

    def __eq__(self, other) -> bool:
        if not isinstance(other, LocalAlignment):
            return False

        return (self.a == other.a and self.b == other.b and
                self.strand == other.strand and self.arange == other.arange
                and self.brange == other.brange)

    def __hash__(self) -> int:
        return hash((self.a, self.b, self.strand, self.arange, self.brange))

    def __len__(self) -> int:
        return self.get_overlap_length()

    def __repr__(self) -> str:
        return ("LocalAlignment[a={0.a_id}, b={0.b_id}, arange={0.arange}, "
                "brange={0.brange}, strand={0.strand}]".format(self))
