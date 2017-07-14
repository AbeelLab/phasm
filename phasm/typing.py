"""
Several type aliases used throughout PHASM
------------------------------------------

This is a separate module to prevent circular imports.
"""

from typing import (Mapping, Set, Callable, Union, Tuple, Iterable, NamedTuple,
                    Any)

# Pairwise local alignments
OrientedDNASegment = 'phasm.alignments.OrientedDNASegment'
OrientedRead = 'phasm.alignments.OrientedRead'
LocalAlignment = 'phasm.alignments.LocalAlignment'
AlignmentsT = Mapping[OrientedRead, Set[LocalAlignment]]

# Assembly Graphs
AssemblyGraph = 'phasm.assembly_graph.AssemblyGraph'
Node = Union[OrientedDNASegment, str]
Edge = Union[Tuple[Node, Node], Tuple[Node, Node, Any]]
Path = Iterable[Edge]
Bubble = Tuple[Node, Node]

# Phasing algorithm parameters
PruneParam = Union[float, Callable[[float], float]]
RelevantReadInfo = NamedTuple('RelevantReadInfo', [
    ('alignments', Set[LocalAlignment]), ('overlap_max', int)
])
RelevantReads = Mapping[OrientedRead, RelevantReadInfo]
