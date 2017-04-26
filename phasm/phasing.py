"""
Haplotype phasing using an assembly graph
=========================================

This module contains the algorithm to phase a haplograph to a set of linear DNA
strings, each representing a chromosome.
"""

import logging
import itertools
from typing import Union, Callable, List, Iterable
from collections import OrderedDict, deque

import networkx

from phasm.alignments import Read
from phasm.assembly_graph import AssemblyGraph
from phasm.bubbles import find_superbubbles, Node

PruneParam = Union[float, Callable[[int], float]]

logger = logging.getLogger(__name__)


class PhasingError(Exception):
    pass


class HaplotypeSet:
    """Represents a set of (incomplete) haplotypes, or a set of `k` DNA strings
    (spelled by the corresponding nodes) each representing a chromosome
    copy."""

    def __init__(self, ploidy: int):
        self.haplotypes = [deque()] * ploidy

    def extend(self, extension: List[List[Read]]):
        pass


class HaploGraphPhaser:
    def __init__(self, g: AssemblyGraph, ploidy: int, threshold: PruneParam,
                 prune_factor: PruneParam):
        self.g = g
        self.ploidy = ploidy
        self.threshold = threshold
        self.prune_factor = prune_factor

        # We start with a single empty haplotype set
        self.possible_haplotypes = [HaplotypeSet(ploidy)]

    def phase(self):
        # find_superbubbles finds superbubbles in reverse topological order
        bubbles = OrderedDict(b for b in find_superbubbles(self.g))
        bubble_sources = set(bubbles.keys())

        start_points = [n for n in self.g if self.g.in_degree(n) == 0]

        if len(start_points) > 1:
            raise PhasingError("There are multiple nodes with in-degree 0, "
                               "not sure what the start point is.")

        curr_node = start_points[0]
        while curr_node:
            if curr_node in bubble_sources:
                # Update possible haplotype sets
                self.branch(curr_node, bubbles[curr_node])
                self.prune(self.prune_factor)

                # Move to the bubble sink
                curr_node = bubbles[curr_node]
            else:
                # This should be a node with a single successor
                if self.g.out_degree(curr_node) > 1:
                    raise PhasingError("Encountered a node which is not a "
                                       "bubble source and has out-degree > 1")
                elif self.g.out_degree(curr_node) == 1:
                    # TODO: extend haplotype sets with single node
                    curr_node = self.g.neighbors(curr_node)[0]
                else:
                    # Last node, no neighbours, quit
                    curr_node = None

        # Only keep the most probable haplotype set
        self.prune(1.0)
        return self.possible_haplotypes

    def branch(self, source: Node, sink: Node):
        """Generate all new possible haplotype extensions for a
        bubble (`source`, `sink`)."""

        new_haplotype_sets = []
        for haplotype_set in self.possible_haplotypes:
            new_haplotype_sets.extend(self.extend_bubble(haplotype_set, source,
                                                         sink))

        self.possible_haplotypes = new_haplotype_sets

    def extend_bubble(self, haplotype_set: HaplotypeSet, source: Node,
                      sink: Node) -> Iterable[HaplotypeSet]:
        """For a given haplotype set, generate all possible extensions at the
        given bubble."""

        possible_paths = networkx.all_simple_paths(self.g, source, sink)

        for extension in itertools.product(possible_paths, repeat=self.ploidy):
            haplotype_set.extend(extension)
            haplotype_set.calculate_rl()

            if haplotype_set.relative_likelihood > self.threshold:
                yield haplotype_set

    def prune(self, prune_factor: PruneParam):
        """Prune the list of possible haplotypes, by removing unlikely
        haplotypes, compared to the most likely haplotype."""

        if len(self.possible_haplotypes) < 2:
            return

        max_prob = max(hs.relative_likelihood for hs in
                       self.possible_haplotypes)

        threshold = (prune_factor(len(self.possible_haplotypes))*max_prob if
                     callable(prune_factor) else prune_factor*max_prob)

        self.possible_haplotypes = [hs for hs in self.possible_haplotypes if
                                    hs.relative_likelihood > threshold]
