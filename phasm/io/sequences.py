"""
Sequence I/O: retreive the actual sequence for a read
=====================================================

This module provides classes and functions to get a the actual sequence
of a read. Currently it's only possible to read them from a FASTA file,
but in the future it could be expanded to obtain these sequences from a
database, and possibly with smarter memory management using an LRU cache.
"""

from abc import ABCMeta, abstractmethod
from itertools import zip_longest

import dinopy

from phasm.alignments import OrientedRead, OrientedDNASegment, MergedReads


class SequenceSource(metaclass=ABCMeta):
    """Base class for any source of sequence data."""

    def __init__(self):
        self.cache = {}  # type: Mapping[OrientedDNASegment, bytes]

    @abstractmethod
    def _get_sequence(self, read: OrientedRead) -> bytes:
        """Get the sequence of the given read, and return as a bytes object.
        This function should take the reverse complement if necessary.

        .. seealso:: get_sequence"""

    def get_sequence(self, read: OrientedDNASegment) -> bytes:
        """Retreive the sequence for a given `OrientedDNASegment`.

        This function also takes care of `MergedReads`, taking the reverse
        complement if necessary, and caches all reads, to make sure we don't
        need to read from disk each time we need a sequence."""

        if read not in self.cache:
            if isinstance(read, MergedReads):
                seq = self._get_merged_reads_sequence(read)
            else:
                seq = self._get_sequence(read)

            self.cache[read] = seq

        return self.cache[read]

    def _get_merged_reads_sequence(self, merged_reads: MergedReads):
        if len(merged_reads.reads) != len(merged_reads.prefix_lengths)+1:
            raise ValueError(
                "Invalid `MergedReads` data structure, not enough information "
                "on read prefix lengths."
            )

        sequence_parts = []
        for read, prefix_len in zip_longest(merged_reads.reads,
                                            merged_reads.prefix_lengths):
            sequence = self.get_sequence(read)
            if prefix_len:
                sequence_parts.append(sequence[:prefix_len])
            else:
                sequence_parts.append(sequence)

        return b"".join(sequence_parts)


class FastaSource(SequenceSource):
    """Obtain sequences from a FASTA file. Creates a FASTA index file for
    faster random access.

    Keeps all read sequences in memory for quick retreival the next time."""

    def __init__(self, file_source):
        super().__init__()
        self.reader = dinopy.FastaReader(file_source)

    def _get_sequence(self, read: OrientedRead) -> bytes:
        """Read the sequence from the FASTA file. Returns the reverse
        complement if necessary."""

        seq = list(self.reader[read.id])[0].sequence.upper()

        if read.orientation == "-":
            seq = dinopy.reverse_complement(seq)

        return seq
