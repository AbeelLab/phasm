import enum
from typing import List, Tuple, NamedTuple


class Strand(enum.IntEnum):
    SAME = 0
    OPPOSITE = 1


class AlignmentType(enum.IntEnum):
    OVERLAP_AB = 0
    OVERLAP_BA = 1
    A_CONTAINED = 2
    B_CONTAINED = 3


_Read = NamedTuple(
    'Read',
    [('id', str), ('moviename', str), ('well', str), ('pulse_start', int),
     ('pulse_end', int), ('length', int)]
)


class Read(_Read):
    def __len__(self) -> int:
        return self.length

    def __str__(self) -> str:
        return self.id

    def __hash__(self) -> int:
        return hash(self.id)

    def get_oriented_read(self, orientation):
        return OrientedRead(self, orientation)


_OrientedRead = NamedTuple('OrientedRead',
                           [('read', Read), ('orientation', str)])


class OrientedRead(_OrientedRead):
    def __getattr__(self, key):
        data = self.read._asdict()

        if key in data:
            return data[key]
        else:
            raise AttributeError("OrientedRead has no attribute '{}'".format(
                key))

    def __str__(self):
        return str(self.read) + self.orientation

    def __hash__(self):
        return hash(str(self.read) + self.orientation)


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
        return self.a.get_oriented_read('+'), self.b.get_oriented_read(
            '+' if self.strand == Strand.SAME else '-')

    def __len__(self):
        return self.get_overlap_length()
