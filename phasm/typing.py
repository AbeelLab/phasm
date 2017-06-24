"""
Several type aliases used througout PHASM
-----------------------------------------

This is a separate module to prevent circular imports.
"""

from typing import Mapping, Set, Callable, Union, Tuple, Iterable

# Pairwise local alignments
OrientedDNASegment = 'phasm.alignments.OrientedDNASegment'
OrientedRead = 'phasm.alignments.OrientedRead'
LocalAlignment = 'phasm.alignments.LocalAlignment'
AlignmentsT = Mapping[OrientedRead, Set[LocalAlignment]]

# Assembly Graphs
AssemblyGraph = 'phasm.assembly_graph.AssemblyGraph'
Node = Union[OrientedDNASegment, str]
Edge = Tuple[Node, Node]
Path = Iterable[Edge]
Bubble = Tuple[Node, Node]

# Phasing algorithm parameters
PruneParam = Union[float, Callable[[float], float]]
RelevantLA = Mapping[OrientedRead, Set[LocalAlignment]]
