import logging
from abc import ABCMeta, abstractmethod

from phasm.alignments import LocalAlignment, AlignmentType


class AlignmentFilter(metaclass=ABCMeta):
    def __init__(self):
        self._filtered = 0
        self.logger = logging.getLogger("{}.{}".format(
            __name__, self.__class__.__name__))

    def __call__(self, la):
        result = self.filter(la)
        if not result:
            self.logger.debug("Filtered local alignment: %s", la)
            self._filtered += 1

        return result

    @property
    def filtered(self):
        """Returns the number of alignments removed."""
        return self._filtered

    @abstractmethod
    def filter(self, la: LocalAlignment):
        pass


class MinReadLength(AlignmentFilter):
    """Factory to generate a filter function which ignores all local alignments
    if one of the reads is shorter than the given minimum read length."""

    def __init__(self, min_read_length: int):
        super().__init__()
        self.min_read_length = min_read_length

    def filter(self, la: LocalAlignment):
        return (len(la.a) >= self.min_read_length and
                len(la.b) >= self.min_read_length)


class MinOverlapLength(AlignmentFilter):
    """Factory to generate a filter function which ignores all local alignments
    with a smaller overlap length than `min_length`.

    To be used with Python's builtin `filter` function, and some iterable which
    outputs `LocalAlignment`s.

    .. seealso:: LocalAlignment, MaxDifferences, MaxErrorRate"""

    def __init__(self, min_length: int):
        super().__init__()
        self.min_length = min_length

    def filter(self, la: LocalAlignment):
        return la.get_overlap_length() >= self.min_length


class ContainedReads(AlignmentFilter):
    """This filter removes all alignments for which holds that one read is
    contained by the other.

    To be used with Python's builtin `filter` function, and some iterable which
    outputs `LocalAlignment`s.

    .. seealso:: LocalAlignment
    """

    def __init__(self):
        super().__init__()
        self.contained_reads = set()

    def filter(self, la: LocalAlignment):
        la_type = la.classify()

        if la_type == AlignmentType.A_CONTAINED:
            self.contained_reads.add(la.a)
            return False
        elif la_type == AlignmentType.B_CONTAINED:
            self.contained_reads.add(la.b)
            return False
        elif la.a in self.contained_reads or la.b in self.contained_reads:
            return False

        return True


class MaxOverhang(AlignmentFilter):
    """This filter removes alignments with too much overhang.

    To be used with Python's builtin `filter` function, and some iterable which
    outputs `LocalAlignment`s.

    .. seealso:: LocalAlignment
    """

    def __init__(self, max_overhang: int, overlap_overhang_ratio: float):
        super().__init__()
        self.max_overhang = max_overhang
        self.ratio = overlap_overhang_ratio

    def filter(self, la: LocalAlignment):
        threshold = min(self.max_overhang, self.ratio*la.get_overlap_length())

        return la.get_overhang() <= threshold
