"""
Haplotype phasing using an assembly graph
=========================================

This module contains the algorithm to phase a haplograph to a set of linear DNA
strings, each representing a chromosome copy.
"""

import math
import random
import logging
import itertools
from typing import Union, Callable, List, Iterable, Mapping, Set
from collections import OrderedDict, deque, defaultdict

import numpy
import networkx

from phasm.alignments import (OrientedRead, LocalAlignment, MergedReads,
                              AlignmentType)
from phasm.assembly_graph import AssemblyGraph
from phasm.bubbles import find_superbubbles, Node, superbubble_nodes


logger = logging.getLogger(__name__)

AlignmentsT = Mapping[OrientedRead, Set[LocalAlignment]]
ReadPositionT = Mapping[OrientedRead, int]
AnchorPointsT = Mapping[int, LocalAlignment]
RelevantAlignmentsT = Mapping[OrientedRead, AnchorPointsT]

ErrorModel = Callable[[numpy.array, numpy.array], float]
PruneParam = Union[float, Callable[[float], float]]


class PhasingError(Exception):
    pass


class MismatchErrorModel(ErrorModel):
    def __init__(self, error_prob: float):
        self.error_prob = error_prob

    def __call__(self, haplo_seq: numpy.array, read_seq: numpy.array) -> float:
        matches = haplo_seq == read_seq
        num_matches = numpy.sum(matches)
        num_mismatches = len(read_seq) - num_matches

        return ((1-self.error_prob)**num_matches *
                self.error_prob**num_mismatches)


class HaplotypeSet:
    """Represents a set of (incomplete) haplotypes, or a set of `k` DNA strings
    (spelled by the corresponding nodes) each representing a chromosome
    copy."""

    def __init__(self, ploidy: int):
        # The nodes for each haplotype
        self.haplotypes = [deque()] * ploidy

        # The actual DNA sequence of each haplotype
        self.haplotypes_dna = [bytearray()] * ploidy

        # Store for the reads spelling the haplotype DNA sequence where its
        # contribution to the DNA sequence starts
        self.read_contribs = [{}] * ploidy  # type: List[ReadPositionT]

        # Log likelihood of this haplotype set
        self.log_likelihood = 0.0

        # For each haplotype keep a list of relevant reads (local alignments)
        # which will be used for likelihood calculation
        self.relevant_reads = defaultdict(dict)  # type: RelevantAlignmentsT

    def extend(self, g: AssemblyGraph, alignments: AlignmentsT,
               extensions: List[List[Node]], last_bubble: bool=False):
        """Extend the haplotype set with a new set of paths."""

        for hap_num, (haplotype_nodes, extension) in enumerate(
                zip(self.haplotypes, extensions)):

            # Determine which read is responsible for which DNA segment of our
            # extension
            pos = len(self.haplotypes_dna[hap_num])
            self.read_contribs[hap_num].update(
                self._get_read_contribs(g, pos, extension,
                                        include_last=last_bubble)
            )

            # Spell the actual DNA sequence
            new_dna = g.sequence_for_path(
                g.node_path_edges(extension, data=True),
                edge_len=g.edge_len,
                include_last=last_bubble
            )

            self.haplotypes_dna[hap_num].extend(new_dna)

            # Add the nodes of the extension to each haplotype
            if (len(haplotype_nodes) > 0 and
                    haplotype_nodes[-1] == extension[0]):
                haplotype_nodes.extend(extension[1:])
            else:
                haplotype_nodes.extend(extension)

            # Update the set of relevant reads, later to be used for likelihood
            # calculation
            for node in extension:
                if isinstance(node, MergedReads):
                    for read in node.reads:
                        self._add_alignments_of_read(alignments, read,
                                                     hap_num)
                else:
                    self._add_alignments_of_read(alignments, node,
                                                 hap_num)

    def _get_read_contribs(self, g: AssemblyGraph, start_pos: int,
                           path: List[Node],
                           include_last: bool=True) -> ReadPositionT:
        """Determine which read is responsible for which piece of DNA in
        the given haplotype"""
        read_start_pos = {}
        total = start_pos
        for u, v, l in g.node_path_edges(path, data=g.edge_len):
            if isinstance(u, MergedReads):
                for read, prefix_len in zip(u.reads, u.prefix_lengths):
                    read_start_pos[read] = total
                    total += prefix_len
            else:
                read_start_pos[u] = total
                total += prefix_len

        if include_last:
            if isinstance(path[-1], MergedReads):
                for read, prefix_len in zip(path[-1].reads,
                                            path[-1].prefix_lengths):
                    read_start_pos[read] = total
                    total += prefix_len
            else:
                read_start_pos[path[-1]] = total

        return read_start_pos

    def _add_alignments_of_read(self, alignments: AlignmentsT,
                                read: OrientedRead, haplotype_num: int):
        """Add all pairwise alignments of a read part of our assembly graph
        to the set of relevant reads needed for likelihood calculation.

        For each alignment, we also store the best anchor point on each
        haplotype."""

        if read not in alignments:
            return

        for la in alignments[read]:
            a_read, b_read = la.get_oriented_reads()

            # The a-read is a read that's part of the haplotype graph.
            # The b-read is a read that aligns to this haplotype.
            # We want to keep all relevant reads (also reads that are
            # not part of the graph), so we store the b-read, with the
            # corresponding best overlapping anchor point
            if (b_read not in self.relevant_reads or
                    haplotype_num not in self.relevant_reads[b_read]):
                self.relevant_reads[b_read][haplotype_num] = la
            else:
                # Only change anchor point if it is better overlapping
                existing_la = self.relevant_reads[b_read][haplotype_num]
                if existing_la.get_overlap_length() < la.get_overlap_length():
                    self.relevant_reads[b_read][haplotype_num] = la

    def calculate_rl(self, g: AssemblyGraph, error_model: ErrorModel,
                     estimated_length: int):
        self.log_likelihood = 0.0

        for read, anchor_points in self.relevant_reads.items():
            read_dna = g.get_sequence(read)
            read_prob = self.calculate_read_prob(read, read_dna, anchor_points,
                                                 error_model, estimated_length)

            if read_prob > 0.0:
                self.log_likelihood += math.log10(read_prob)

    def calculate_read_prob(self, read: OrientedRead, read_dna: bytes,
                            anchor_points: AnchorPointsT,
                            error_model: ErrorModel,
                            estimated_length: int) -> float:
        read_prob = 0.0

        for haplotype_num, la in anchor_points.items():
            haplotype_dna = numpy.frombuffer(
                self.haplotypes_dna[haplotype_num], dtype='b')
            read_dna = numpy.frombuffer(read_dna, dtype='b')
            read_contribs = self.read_contribs[haplotype_num]

            a_read, b_read = la.get_oriented_reads()
            anchor_start = read_contribs[a_read]
            hseq_start = self._get_alignment_start_pos(anchor_start, la)

            # If the read is longer and extends further than our current
            # assembled haplotype, then ignore these additional bases for
            # now. These additional bases will be phased later when we move
            # on to the next bubble.
            seq_length = min(len(haplotype_dna[hseq_start:]),
                             len(read_dna))

            # Get the subsequences which we're going to compare
            read_subseq = read_dna[:seq_length]
            haplo_subseq = haplotype_dna[
                hseq_start:hseq_start+seq_length]

            prob = error_model(haplo_subseq, read_subseq)
            prob *= 1/estimated_length

            read_prob += prob

        return read_prob

    def _get_alignment_start_pos(self, anchor_start: int,
                                 alignment: LocalAlignment) -> int:
        la_type = alignment.classify()

        # The a-read should be part of the haplotype graph, and is used to
        # anchor any other alignments to this read.
        if la_type in (AlignmentType.OVERLAP_AB, AlignmentType.B_CONTAINED):
            # Haplotype: --------------------------------------------------
            #    A read:           |------|----->
            #    B read:                  |-------------->
            #                           /    (possibly some overhang)
            return max(
                0, anchor_start + alignment.arange[0] - alignment.brange[0]
            )
        elif la_type in (AlignmentType.OVERLAP_BA, AlignmentType.A_CONTAINED):
            # Haplotype: --------------------------------------------------
            #    A read:           |------------>
            #    B read:    |------|-------->
            return max(0, anchor_start - alignment.brange[0])


class HaploGraphPhaser:
    def __init__(self, g: AssemblyGraph,
                 alignments: Mapping[Node, Set[LocalAlignment]],
                 ploidy: int,
                 error_model: ErrorModel,
                 threshold: PruneParam,
                 prune_factor: PruneParam):
        self.g = g
        self.alignments = alignments
        self.error_model = error_model
        self.ploidy = ploidy
        self.num_bubbles = 0

        self.threshold = threshold
        self.prune_factor = prune_factor

        # We start with a single empty haplotype set
        self.possible_haplotypes = [HaplotypeSet(ploidy)]

    def phase(self):
        # find_superbubbles finds superbubbles in reverse topological order
        bubbles = OrderedDict(b for b in find_superbubbles(self.g))
        bubble_sources = set(bubbles.keys())
        self.num_bubbles = len(bubbles)

        start_points = [n for n in self.g if self.g.in_degree(n) == 0]
        end_points = [n for n in self.g if self.g.out_degree(n) == 0]

        if len(start_points) > 1:
            raise PhasingError("There are multiple nodes with in-degree 0, "
                               "not sure what the start point is.")

        if len(end_points) > 1:
            raise PhasingError("Found multiple 'sink' nodes, this should not "
                               "be the case.")

        # Calculate approximate haplotype lengths, by averaging the length of
        # a few random paths from start to end
        self.approx_length = self._calculate_approx_length(bubbles)

        logger.info("Start phasing haplograph with %d bubbles and approximated"
                    " length of %d bases.", self.num_bubbles,
                    self.approx_length)

        curr_node = start_points[0]
        bubble_num = 0
        while curr_node:
            if curr_node in bubble_sources:
                # Update possible haplotype sets
                logger.debug("Branch and prune for bubble %d/%d",
                             bubble_num+1, self.num_bubbles)
                self.branch(curr_node, bubbles[curr_node], bubble_num)
                self.prune(self.prune_factor, bubble_num)

                # Move to the bubble sink
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

    def branch(self, source: Node, sink: Node, bubble_num: int):
        """Generate all new possible haplotype extensions for a
        bubble (`source`, `sink`)."""

        bubble_alignments = set()
        for node in superbubble_nodes(self.g, source, sink):
            if node in self.alignments:
                bubble_alignments.update(self.alignments[node])

        new_haplotype_sets = []
        for haplotype_set in self.possible_haplotypes:
            new_haplotype_sets.extend(
                self.extend_bubble(haplotype_set, source, sink, bubble_num))

        self.possible_haplotypes = new_haplotype_sets

    def extend_bubble(self, haplotype_set: HaplotypeSet, source: Node,
                      sink: Node, bubble_num: int) -> Iterable[HaplotypeSet]:
        """For a given haplotype set, generate all possible extensions at the
        given bubble."""

        possible_paths = networkx.all_simple_paths(self.g, source, sink)

        threshold = self.threshold
        if callable(threshold):
            threshold = threshold(bubble_num / self.num_bubbles)

        threshold = math.log10(threshold)

        if bubble_num == 0:
            # For the first bubble the order does not matter, as a permutation
            # in a different order will in the end result in the same haplotype
            # set
            extension_iter = iter(
                itertools.combinations_with_replacement(possible_paths,
                                                        self.ploidy)
            )
        else:
            extension_iter = iter(itertools.product(possible_paths,
                                                    repeat=self.ploidy))

        for extension in extension_iter:
            haplotype_set.extend(self.g, self.alignments, extension,
                                 bubble_num == self.num_bubbles-1)
            haplotype_set.calculate_rl(self.g, self.error_model,
                                       self.approx_length)

            if haplotype_set.relative_likelihood > threshold:
                yield haplotype_set

    def prune(self, prune_factor: PruneParam, bubble_num: int):
        """Prune the list of possible haplotypes, by removing unlikely
        haplotypes, compared to the most likely haplotype."""

        if len(self.possible_haplotypes) < 2:
            return

        num_before = len(self.possible_haplotypes)

        if callable(prune_factor):
            prune_factor = prune_factor(bubble_num / (self.num_bubbles-1))

        prune_factor = math.log10(prune_factor)

        max_likelihood = max(hs.log_likelihood for hs in
                             self.possible_haplotypes)

        self.possible_haplotypes = [
            hs for hs in self.possible_haplotypes if
            (hs.log_likelihood - max_likelihood) >= prune_factor
        ]

        num_after = len(self.possible_haplotypes)

        logger.debug("Pruned %d unlikely haplotype sets, %d left.",
                     (num_before - num_after), num_after)

    def _calculate_approx_length(self, bubbles: Mapping[Node, Node],
                                 num_paths: int=3) -> int:
        """Calculate the average length of an haplotype by averaging the length
        of `num_path` random paths from start to finish."""

        lengths = []
        for i in range(num_paths):
            path = []
            logger.debug("Generating random path for graph with %d bubbles.",
                         self.num_bubbles)
            for j, (b_src, b_sink) in enumerate(reversed(bubbles.items())):
                possible_paths = list(networkx.all_simple_paths(self.g, b_src,
                                                                b_sink))
                path_extension = random.choice(possible_paths)

                if len(path) > 0 and path[-1] == path_extension[0]:
                    path.extend(path_extension[1:])
                else:
                    path.extend(path_extension)

            logger.debug("Path: %s", path)
            length = sum(l for u, v, l in
                         self.g.node_path_edges(path, data=self.g.edge_len))
            length += len(path[-1])
            lengths.append(length)

        return int(round(sum(lengths) / num_paths))
