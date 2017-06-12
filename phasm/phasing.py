"""
Haplotype phasing using an assembly graph
=========================================

This module contains the algorithm to phase a haplograph to a set of linear DNA
strings, each representing a chromosome copy.
"""

import random
import logging
from itertools import combinations_with_replacement, product
from typing import Union, Callable, List, Iterable, Mapping, Set, Tuple
from collections import OrderedDict, deque

import numpy
import networkx

from phasm.alignments import (OrientedRead, OrientedDNASegment, LocalAlignment,
                              MergedReads)
from phasm.assembly_graph import AssemblyGraph
from phasm.bubbles import find_superbubbles, superbubble_nodes


logger = logging.getLogger(__name__)

Node = OrientedDNASegment
AlignmentsT = Mapping[OrientedRead, Set[LocalAlignment]]
ReadPositionT = Mapping[OrientedRead, int]
AnchorPointsT = Mapping[int, LocalAlignment]
RelevantAlignmentsT = Mapping[OrientedRead, AnchorPointsT]

ErrorModel = Callable[[numpy.array, numpy.array], float]
PruneParam = Union[float, Callable[[float], float]]


class PhasingError(Exception):
    pass


class HaplotypeSet:
    """Represents a set of (incomplete) haplotypes, or a set of `k` DNA strings
    (spelled by the corresponding nodes) each representing a chromosome
    copy."""

    def __init__(self, ploidy: int, copy_from: 'HaplotypeSet'=None):
        self.ploidy = ploidy

        # Nodes spelling each haplotype
        self.haplotypes = []  # type: List[List[Node]]

        # Also keep a set of reads used for each haplotype, useful for
        # relative likelihood calculation
        self.read_set = []  # type: List[Set[OrientedRead]]

        if isinstance(copy_from, HaplotypeSet):
            for i in range(ploidy):
                self.haplotypes.append(deque(copy_from.haplotypes[i]))
                self.read_set.append(set(copy_from.read_set[i]))
        else:
            for i in range(ploidy):
                self.haplotypes.append(deque())
                self.read_set.append(set())

        self.relative_likelihood = 0.0

    def extend(self, extensions: List[Tuple[Node]],
               ext_read_sets: List[Set[OrientedRead]]) -> 'HaplotypeSet':
        """Extend the haplotype set with a new set of paths."""

        # Make a copy of itself for a new set
        new_set = HaplotypeSet(self.ploidy, copy_from=self)

        for hap_num, (extension, read_set) in enumerate(
                zip(extensions, ext_read_sets)):
            haplotype_nodes = new_set.haplotypes[hap_num]

            # Add the nodes of the extension to each haplotype
            # It's possible that the last node of this haplotype set
            # (which is probably a bubble exit), is also the bubble entrance
            # and thus our start node of our extension.
            if (len(haplotype_nodes) > 0 and
                    haplotype_nodes[-1] == extension[0]):
                haplotype_nodes.extend(extension[1:])
            else:
                haplotype_nodes.extend(extension)

            new_set.read_set[hap_num].update(read_set)

            logger.debug("H%d: %s", hap_num, haplotype_nodes)

        return new_set


class HaploGraphPhaser:
    def __init__(self, g: AssemblyGraph,
                 alignments: AlignmentsT,
                 ploidy: int,
                 threshold: PruneParam,
                 prune_factor: PruneParam):
        self.g = g
        self.alignments = alignments
        self.ploidy = ploidy

        self.num_bubbles = 0
        self.prev_bubble = None

        # Keep a record of all alignments seen
        self.prev_alignments = set()

        self.threshold = threshold
        self.prune_factor = prune_factor

        # We start with a single empty haplotype set
        self.possible_haplotypes = [HaplotypeSet(ploidy)]

    def phase(self):
        # find_superbubbles finds superbubbles in reverse topological order
        bubbles = OrderedDict(b for b in
                              find_superbubbles(self.g, report_nested=False))
        bubble_sources = set(bubbles.keys())
        self.num_bubbles = len(bubbles)

        start_points = [n for n in self.g.nodes_iter()
                        if self.g.in_degree(n) == 0]
        end_points = [n for n in self.g.nodes_iter()
                      if self.g.out_degree(n) == 0]

        if len(start_points) > 1:
            raise PhasingError("There are multiple nodes with in-degree 0, "
                               "not sure what the start point is.")

        if len(end_points) > 1:
            raise PhasingError("Found multiple 'sink' nodes, this should not "
                               "be the case.")

        # Calculate approximate haplotype lengths, by averaging the length of
        # a few random paths from start to end
        logger.debug("Superbubbles: %s", bubbles)

        logger.info("Start phasing haplograph with %d bubbles",
                    self.num_bubbles)

        curr_node = start_points[0]
        bubble_num = 0
        while curr_node:
            if curr_node in bubble_sources:
                # Update possible haplotype sets
                logger.info("Branch and prune for bubble %d/%d",
                            bubble_num+1, self.num_bubbles)
                self.branch(curr_node, bubbles[curr_node], bubble_num)
                self.prune(self.prune_factor, bubble_num)

                # Move to the bubble exit
                self.prev_bubble = (curr_node, bubbles[curr_node])
                curr_node = bubbles[curr_node]
                bubble_num += 1
            else:
                # This should be a node with a single successor
                if self.g.out_degree(curr_node) > 1:
                    raise PhasingError("Encountered a node which is not a "
                                       "bubble source and has out-degree > 1")
                elif self.g.out_degree(curr_node) == 1:
                    logger.waring("Encountered a node that is not part of a"
                                  "bubble.")
                    # TODO: extend haplotype sets with single node
                    curr_node = self.g.neighbors(curr_node)[0]
                else:
                    # Last node, no neighbours, quit
                    curr_node = None

        # Only keep the most probable haplotype set
        self.prune(1.0, self.num_bubbles-1)
        logger.info("Done")
        return self.possible_haplotypes

    def branch(self, entrance: Node, exit: Node, bubble_num: int):
        """Generate all new possible haplotype extensions for a
        bubble (`entrance`, `exit`)."""

        if len(self.possible_haplotypes) == 0:
            logger.warning("No possible haplotypes available.")

        # Determine alignments that span between paths in previous and the
        # current bubble.
        current_bubble_alignments = self.get_interior_alignments(entrance,
                                                                 exit)

        spanning_alignments = (self.prev_alignments
                               & current_bubble_alignments)

        if self.prev_bubble:
            logger.info("%d spanning alignments between bubble <%s, %s> and "
                        "bubble <%s, %s>", len(spanning_alignments),
                        *self.prev_bubble, entrance, exit)

        if len(spanning_alignments) == 0 and len(self.possible_haplotypes) > 1:
            # No alignments spanning the previous bubbles and the current one.
            # This means we're actually starting a new haplotype block. But,
            # our `possible_haplotypes` may still be filled with possible
            # candidates, from which none we can link to the current bubble.
            # So, we only keep the most likely candidates, and pick a random
            # one from that list.
            self.prune(1.0, bubble_num)
            random_candidate = random.choice(self.possible_haplotypes)
            self.possible_haplotypes = [random_candidate]

        # We're using a list of tuples because tuples are immutable and can be
        # used as dict keys, and more importantly in a `collections.Counter` to
        # determine the multiplicity of path in a later stage.
        possible_paths = list(
            tuple(p) for p in networkx.all_simple_paths(self.g, entrance, exit)
        )

        logger.info("%d possible paths through bubble <%r, %r>",
                    len(possible_paths), entrance, exit)

        new_haplotype_sets = []
        for haplotype_set in self.possible_haplotypes:
            new_haplotype_sets.extend(
                self.generate_new_hsets(haplotype_set, possible_paths,
                                        spanning_alignments, bubble_num)
            )

        self.possible_haplotypes = new_haplotype_sets
        self.prev_bubble = (entrance, exit)
        self.prev_alignments.update(current_bubble_alignments)

    def generate_new_hsets(self, haplotype_set: HaplotypeSet,
                           possible_paths: List[Tuple[OrientedDNASegment]],
                           spanning_alignments: Set[OrientedRead],
                           bubble_num: int) -> Iterable[HaplotypeSet]:
        """For a given haplotype set, generate all possible extensions at the
        current bubble. Only yield the new haplotype sets that have a relative
        likelihood above a given threshold."""

        if callable(self.threshold):
            threshold = self.threshold(bubble_num / self.num_bubbles)
        else:
            threshold = self.threshold

        if bubble_num == 0 or len(spanning_alignments) == 0:
            # For the first bubble the order does not matter, as a permutation
            # in a different order will in the end result in the same haplotype
            # set. The same holds when there are no spanning alignments between
            # the previous bubble and the current one. We have no way of
            # linking previous haplotypes to the current one (no read evidence)
            # so this is like starting a new haplotype block
            extension_iter = iter(combinations_with_replacement(possible_paths,
                                                                self.ploidy))
        else:
            # Otherwise all possible k-tuples of possible paths
            extension_iter = iter(product(possible_paths, repeat=self.ploidy))

        for extension in extension_iter:
            ext_read_sets = []
            for hap_ext in extension:
                ext_read_sets.append(self.get_extension_reads(hap_ext))

            if len(spanning_alignments) > 0:
                rl = self.calculate_rl(haplotype_set, extension,
                                       ext_read_sets, spanning_alignments)
            else:
                rl = 1.0

            if rl >= threshold:
                new_set = haplotype_set.extend(extension, ext_read_sets)
                new_set.relative_likelihood = rl
                yield new_set

    def prune(self, prune_factor: PruneParam, bubble_num: int):
        """Prune the list of possible haplotypes, by removing unlikely
        haplotypes, compared to the most likely haplotype."""

        num_before = len(self.possible_haplotypes)
        logger.debug("Number of possible haplotypes: %d", num_before)

        if len(self.possible_haplotypes) < 2:
            return

        if callable(prune_factor):
            prune_factor = prune_factor(bubble_num / (self.num_bubbles-1))

        max_likelihood = max(hs.relative_likelihood for hs in
                             self.possible_haplotypes)

        if prune_factor == 1.0:
            logger.info("Max likelihood: %.3f", max_likelihood)
            logger.info("Other likelihoods: %s", ", ".join(
                str(hs.relative_likelihood)
                for hs in self.possible_haplotypes))

        self.possible_haplotypes = [
            hs for hs in self.possible_haplotypes if
            (hs.relative_likelihood / max_likelihood) >= prune_factor
        ]

        num_after = len(self.possible_haplotypes)

        logger.info("Pruned %d unlikely haplotype sets, %d left.",
                    (num_before - num_after), num_after)

    def get_interior_alignments(self, entrance: Node,
                                exit: Node) -> Set[OrientedRead]:
        """Get all local alignments that involves some node inside the given
        bubble, but ignoring alignments that involve the entrance or exit."""

        entrance_reads = set()
        if isinstance(entrance, MergedReads):
            for read in entrance.reads:
                entrance_reads.add(read)
        else:
            entrance_reads.add(entrance)

        exit_reads = set()
        if isinstance(exit, MergedReads):
            for read in exit.reads:
                exit_reads.add(read)

        alignments = set()

        for node in superbubble_nodes(self.g, entrance, exit):
            if isinstance(node, MergedReads):
                for read in node.reads:
                    alignments.update(self._get_alignments_of_read(
                        entrance_reads, exit_reads, read))
            else:
                alignments.update(self._get_alignments_of_read(
                    entrance_reads, exit_reads, node))

        return alignments

    def _get_alignments_of_read(self, entrance_reads: Set[OrientedRead],
                                exit_reads: Set[OrientedRead],
                                read: OrientedRead) -> Iterable[
                                    LocalAlignment]:
        if read not in self.alignments:
            return

        for la in self.alignments[read]:
            a_read, b_read = la.get_oriented_reads()

            if a_read in entrance_reads or a_read in exit_reads:
                continue

            if b_read in entrance_reads or b_read in exit_reads:
                continue

            yield la

    def get_extension_reads(self, extension: Tuple[Node]) -> Set[OrientedRead]:
        reads = set()
        for node in extension:
            if isinstance(node, MergedReads):
                for read in node.reads:
                    reads.add(read)
            else:
                reads.add(node)

        return reads

    def calculate_rl(self, hs: HaplotypeSet, extension: List[Tuple[Node]],
                     ext_read_sets: List[Set[OrientedRead]],
                     spanning_alignments: Set[LocalAlignment]) -> float:
        """Calculate the relative likelihood of a haplotype set extension,
        assuming the given haplotype set is correct."""

        hs_spanning_alignments = set()
        for la in spanning_alignments:
            a_read, b_read = la.get_oriented_reads()

            for hap_read_set, ext_read_set in zip(hs.read_set, ext_read_sets):
                if a_read not in hap_read_set and a_read not in ext_read_set:
                    continue

                if b_read not in hap_read_set and b_read not in ext_read_set:
                    continue

                hs_spanning_alignments.add(la)
                break

        relative_likelihood = (len(hs_spanning_alignments) /
                               len(spanning_alignments))

        # TODO: Include coverage analysis

        return relative_likelihood
