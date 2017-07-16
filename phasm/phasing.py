"""
Haplotype phasing using an assembly graph
=========================================

This module contains the algorithm to phase a haplograph to a set of linear DNA
strings, each representing a chromosome copy.
"""

import sys
import math
import random
import logging
from itertools import combinations_with_replacement, product
from typing import List, Iterable, Set, Tuple
from collections import OrderedDict, deque, Counter

import networkx
from scipy.stats import norm
from scipy.special import binom

from phasm.alignments import OrientedRead, OrientedDNASegment, MergedReads
from phasm.assembly_graph import AssemblyGraph
from phasm.bubbles import find_superbubbles, superbubble_nodes
from phasm.typing import (Node, AlignmentsT, PruneParam, RelevantReadInfo,
                          RelevantReads)

logger = logging.getLogger(__name__)


class PhasingError(Exception):
    pass


class CoverageModel:
    """When we use a certain path through a bubble multiple times (for multiple
    haplotypes), then you would expect that the average coverage for this path
    is higher.

    This model tries to calculate the probability that a certain path occurs
    in multiple haplotypes by assuming that the average coverage is normally
    distributed, and looks at the difference between the expected coverage and
    the actual average coverage."""

    def __init__(self, mu: float, sigma: float):
        self.mu = mu
        self.sigma = sigma

    def calculate_prob(self, multiplicity: int, avg_coverage: float):
        mu = self.mu * multiplicity
        sigma = self.sigma * multiplicity

        difference = math.fabs(mu - avg_coverage)

        return 2*norm.cdf(mu - difference, loc=mu, scale=sigma)


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
        self.read_sets = []  # type: List[Set[OrientedRead]]

        if isinstance(copy_from, HaplotypeSet):
            for i in range(ploidy):
                self.haplotypes.append(deque(copy_from.haplotypes[i]))
                self.read_sets.append(set(copy_from.read_sets[i]))
        else:
            for i in range(ploidy):
                self.haplotypes.append(deque())
                self.read_sets.append(set())

        self.log_rl = float('-inf')

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

            new_set.read_sets[hap_num].update(read_set)

        return new_set


class BubbleChainPhaser:
    def __init__(self, g: AssemblyGraph,
                 ploidy: int,
                 min_spanning_reads: int,
                 threshold: PruneParam,
                 prune_factor: PruneParam):
        self.g = g
        self.ploidy = ploidy

        self.start_of_block = True

        # Keep a record of all reads seen
        self.prev_aligning_reads = set()  # type: Set[OrientedRead]

        # Keep a record of reads playing a role in the assembly graph
        self.prev_graph_reads = set()  # type: Set[OrientedRead]

        self.min_spanning_reads = min_spanning_reads
        self.threshold = threshold
        self.prune_factor = prune_factor

        self.candidate_sets = [HaplotypeSet(self.ploidy)]

        # find_superbubbles finds superbubbles in reverse topological order
        # Superbubbles without an interior are filtered, to make it a bit
        # easier to identify subgraphs which are just linear paths instead of
        # an actual bubble chain.
        def non_branching_bubble(b):
            """Filter superbubbles with an empty interior (i.e. there are no
            other nodes between the entrance and exit)."""

            s, t = b
            return not (g.out_degree(s) == 1 and g.in_degree(t) == 1 and
                        g.successors(s)[0] == t)

        self.bubbles = OrderedDict(b for b in filter(
            non_branching_bubble,
            find_superbubbles(self.g, report_nested=False)
        ))

        self.bubble_entrances = set(self.bubbles.keys())
        self.num_bubbles = len(self.bubbles)

        self.start_points = [n for n in self.g.nodes_iter()
                             if self.g.in_degree(n) == 0]
        self.end_points = [n for n in self.g.nodes_iter()
                           if self.g.out_degree(n) == 0]

        if len(self.start_points) > 1:
            raise PhasingError("There are multiple nodes with in-degree 0, "
                               "not sure what the start point is.")

        if len(self.end_points) > 1:
            raise PhasingError("Found multiple 'sink' nodes, this should not "
                               "be the case.")

        logger.debug("Superbubbles: %s", self.bubbles)

    def phase(self, alignments: AlignmentsT) -> Iterable[
            Tuple[HaplotypeSet, bool]]:
        logger.info("Start phasing haplograph with %d bubbles",
                    self.num_bubbles)
        self.alignments = alignments

        entrance = self.start_points[0]

        if entrance not in self.bubble_entrances:
            raise PhasingError("Unexpected starting point '{}', this is not "
                               "a bubble entrance".format(entrance))

        bubble_num = 0
        num_bubbles = len(self.bubbles)
        while entrance in self.bubble_entrances:
            exit = self.bubbles[entrance]

            logger.info("-----")
            logger.info("Current bubble <%s, %s> (%d/%d)", entrance, exit,
                        bubble_num+1, num_bubbles)

            # Determine reads that span between paths in previous bubbles and
            # the current bubble.

            # The difference between `cur_bubble_graph_reads` and
            # `cur_bubble_aligning_reads` is that the former one only contains
            # reads that play a role in the assembly graph, and the latter one
            # also contains reads that don't necessarily play a role in the
            # assembly graph, but do align to the current bubble.
            #
            # Additionaly, both variables only include reads that are related
            # to the *interior* of the current bubble, so an alignment to the
            # entrance or exit is ignored. Alignments to the entrance or exit
            # do not provide any information on connecting paths between two
            # different bubbles.
            interior_nodes = (superbubble_nodes(
                self.g, entrance, exit) - {entrance, exit})
            cur_bubble_graph_reads = set(self.get_all_reads(interior_nodes))
            cur_bubble_aligning_reads = self.get_aligning_reads(interior_nodes)
            logger.debug("Current bubble has %d reads used for its interior "
                         "and %d reads aligning to this bubble.",
                         len(cur_bubble_graph_reads),
                         len(cur_bubble_aligning_reads))
            spanning_reads = (self.prev_aligning_reads &
                              cur_bubble_aligning_reads)

            if (len(spanning_reads) < self.min_spanning_reads and not
                    self.start_of_block):
                # Not enough reads to link the current bubble to the previous
                # bubbles, yield current haplotype block and start a new one.
                logger.info("Starting new haplotype block because we have too "
                            "few spanning reads (%d < %d)",
                            len(spanning_reads),
                            self.min_spanning_reads)
                yield self.new_block(), False

            # Based on the relevant reads, collect all relevant local
            # alignments. We ignore local alignments downstream of the
            # current bubble.
            relevant_reads = (spanning_reads if not self.start_of_block else
                              cur_bubble_aligning_reads)
            rel_read_info = self.get_relevant_read_info(
                relevant_reads, cur_bubble_graph_reads)

            total_relevant_la = 0
            for read_alignments, _ in rel_read_info.values():
                total_relevant_la += len(read_alignments)

            logger.info("%d relevant reads inducing %d relevant local "
                        "alignments", len(relevant_reads),
                        total_relevant_la)

            new_haplotype_sets = self.branch(entrance, exit, rel_read_info,
                                             not self.start_of_block)

            if new_haplotype_sets:
                if not self.start_of_block:
                    new_haplotype_sets = self.prune(new_haplotype_sets,
                                                    self.prune_factor)
                self.candidate_sets = new_haplotype_sets
            else:
                # None of the branched haplotype sets are above threshold,
                # yield current haploblock and start a new block at the current
                # bubble.
                logger.info("None of the branched candidates at this bubble "
                            "has a relative likelihood above the given "
                            "threshold. Starting a new block. If this happens "
                            "too often consider lowering your threshold.")
                yield self.new_block(), False

                # Try again as start of a new block (with empty relevant_la).
                # This list is guaranteed to be filled because candidates
                # generated by the first bubble in a block are not yet
                # filtered.
                self.candidate_sets = self.branch(entrance, exit,
                                                  rel_read_info,
                                                  check_threshold=False)

            # Move to next bubble
            entrance = exit
            self.prev_graph_reads.update(cur_bubble_graph_reads)
            self.prev_aligning_reads.update(cur_bubble_aligning_reads)
            bubble_num += 1
            self.start_of_block = False

        # Only keep the most probable haplotype set(s)
        self.candidate_sets = self.prune(self.candidate_sets, 1.0)
        logger.info("Got %d equally likely candidate sets, picking a random "
                    "one.", len(self.candidate_sets))
        yield random.choice(self.candidate_sets), True
        logger.info("Done")

    def new_block(self):
        # Only keep the most probable candidate
        best_candidates = self.prune(self.candidate_sets, 1.0)
        logger.info("Got %d equally likely candidate sets, picking a random "
                    "one.", len(best_candidates))
        random_best = random.choice(best_candidates)

        # Reset state
        self.candidate_sets = [HaplotypeSet(self.ploidy)]
        self.prev_graph_reads = set()
        self.prev_aligning_reads = set()
        self.start_of_block = True

        return random_best

    def get_relevant_read_info(self, relevant_reads: Set[OrientedRead],
                               cur_bubble_graph_reads: Set[OrientedRead]
                               ) -> RelevantReads:
        """Get all alignments of the given spanning reads that either involve
        one of the reads in the interior of one of the previous bubbles or
        one of the reads in the interior of the current bubble.

        This ensures we only obtain alignments that do not involve an entrance
        or exit, and we ignore any alignments in following bubbles.
        """

        relevant_reads_info = {}
        for read in relevant_reads:
            alignments = set()
            overlap_max = 0
            for la in self.alignments[read]:
                a_read, b_read = la.get_oriented_reads()

                if (b_read not in cur_bubble_graph_reads and
                        b_read not in self.prev_graph_reads):
                    # Alignment is part of the entrance or exit,
                    # or the alignment does not involve any reads part of the
                    # assembly graph, or the alignment involves a read further
                    # down the graph (beyond the current bubble)
                    continue

                alignments.add(la)
                overlap_max = max(overlap_max, la.arange[1])

            if alignments:
                logger.debug("Relevant read %s with overlap_max = %d, "
                             "alignments: %s", overlap_max, alignments)
                relevant_reads_info[read] = RelevantReadInfo(alignments,
                                                             overlap_max)

        return relevant_reads_info

    def branch(self, entrance: Node, exit: Node,
               relevant_reads: RelevantReads,
               check_threshold: bool=True) -> List[HaplotypeSet]:
        """Generate all new possible haplotype extensions for a
        bubble (`entrance`, `exit`)."""

        # We're using a list of tuples because tuples are immutable and can be
        # used as dict keys, and more importantly in a `collections.Counter` to
        # determine the multiplicity of path in a later stage.
        possible_paths = list(
            tuple(p) for p in networkx.all_simple_paths(self.g, entrance, exit)
        )

        logger.info("%d possible paths through bubble <%r, %r>",
                    len(possible_paths), entrance, exit)

        new_haplotype_sets = []
        for haplotype_set in self.candidate_sets:
            new_haplotype_sets.extend(
                self.generate_new_hsets(haplotype_set, possible_paths,
                                        relevant_reads, check_threshold)
            )

        return new_haplotype_sets

    def generate_new_hsets(self, haplotype_set: HaplotypeSet,
                           possible_paths: List[Tuple[OrientedDNASegment]],
                           relevant_reads: RelevantReads,
                           check_threshold: bool=True
                           ) -> Iterable[HaplotypeSet]:
        """For a given haplotype set, generate all possible extensions at the
        current bubble. Only yield the new haplotype sets that have a relative
        likelihood above a given threshold."""

        if callable(self.threshold):
            threshold = self.threshold(len(relevant_reads))
        else:
            threshold = self.threshold

        if threshold == 0.0:
            threshold = float('-inf')
        else:
            threshold = math.log10(threshold)

        if self.start_of_block:
            # For the first bubble the order does not matter, as a permutation
            # in a different order will in the end result in the same haplotype
            # set.
            extension_iter = iter(combinations_with_replacement(possible_paths,
                                                                self.ploidy))
            num_possible_sets = binom(self.ploidy + len(possible_paths) - 1,
                                      self.ploidy)
        else:
            # Otherwise all possible k-tuples of possible paths, because now
            # order does matter
            extension_iter = iter(product(possible_paths, repeat=self.ploidy))
            num_possible_sets = self.ploidy**len(possible_paths)

        for extension in extension_iter:
            ext_read_sets = []
            for hap_ext in extension:
                # We index with [1:-1] to ignore the entrance and exit of the
                # bubble
                ext_read_sets.append(set(self.get_all_reads(hap_ext[1:-1])))

            # For the first bubble in a new haploblock we just generate all
            # candidates without any relative likelihood calculation, because
            # we don't have enough spanning reads to do so.
            rl = self.calculate_rl(haplotype_set, extension, ext_read_sets,
                                   relevant_reads, num_possible_sets)

            if not check_threshold:
                new_set = haplotype_set.extend(extension, ext_read_sets)
                new_set.log_rl = rl
                yield new_set
            elif rl >= threshold:
                new_set = haplotype_set.extend(extension, ext_read_sets)
                new_set.log_rl = rl
                yield new_set

    def prune(self, candidate_sets: List[HaplotypeSet],
              prune_factor: PruneParam):
        """Prune the list of candidate haplotype sets, by removing unlikely
        haplotype sets, compared to the most likely haplotype set."""

        num_before = len(candidate_sets)
        logger.debug("Number of candidate sets: %d", num_before)

        if len(candidate_sets) < 2:
            return candidate_sets

        if callable(prune_factor):
            prune_factor = prune_factor(len(candidate_sets))

        if prune_factor <= 0.0:
            return candidate_sets

        prune_factor = math.log10(prune_factor)
        max_likelihood = max(hs.log_rl for hs in candidate_sets)

        logger.info("Max log relative likelihood: %.5f", max_likelihood)
        logger.debug("Other log relative likelihoods: %s", ", ".join(
            "{:.2f}".format(hs.log_rl) for hs in candidate_sets))

        if max_likelihood != float('-inf'):
            candidate_sets = [
                hs for hs in candidate_sets if
                (hs.log_rl - max_likelihood) >= prune_factor
            ]

        num_after = len(candidate_sets)

        logger.info("Pruned %d unlikely haplotype sets with prune factor %.2f,"
                    " %d left.",
                    (num_before - num_after), 10.0**prune_factor, num_after)

        return candidate_sets

    def get_aligning_reads(self, nodes: List[Node]) -> Set[OrientedRead]:
        aligning_reads = set()  # type: Set[OrientedRead]

        for read in self.get_all_reads(nodes):
            aligning_reads.add(read)
            if read in self.alignments:
                for la in self.alignments[read]:
                    a_read, b_read = la.get_oriented_reads()
                    aligning_reads.add(b_read)

        return aligning_reads

    def get_all_reads(self, nodes: Iterable[Node]) -> Iterable[OrientedRead]:
        """Collects all relevant `OrientedRead` objects for a given list of
        nodes from the assembly graph. Basically it expands any encountered
        `MergedReads` instance."""

        for node in nodes:
            if isinstance(node, MergedReads):
                for read in node.reads:
                    yield read
            else:
                yield node

    def calculate_rl(self, hs: HaplotypeSet, extension: List[Tuple[Node]],
                     ext_read_sets: List[Set[OrientedRead]],
                     relevant_reads: RelevantReads,
                     num_possible_sets: int) -> float:
        """Calculate the relative likelihood of a haplotype set extension,
        assuming the given haplotype set is correct."""

        # TODO: Coverage model is not good enough yet, disable for now
        #
        # Count how often a path is used in this haplotype set
        # multiplicities = Counter(extension)

        # coverage_prob = 1.0
        # include_last = bubble_num == self.num_bubbles - 1
        # for path in extension:
        #     avg_coverage = average_coverage_path(
        #         self.g, self.alignments, path, include_last)

        #     coverage_prob *= self.coverage_model.calculate_prob(
        #         multiplicities[path], avg_coverage)

        # logger.debug("Coverage probability: %.3f", coverage_prob)
        # if len(relevant_la) == 0:
        #     return coverage_prob

        if len(relevant_reads) < self.min_spanning_reads:
            return 0.0

        logger.debug("Haplotype read sets: %s", hs.read_sets)
        logger.debug("Extension read sets: %s", ext_read_sets)

        # Calculate P[R|H,E]
        ext_prob = 0.0
        for read, info in relevant_reads.items():
            read_prob = 0.0
            for hap_read_set, ext_read_set in zip(hs.read_sets, ext_read_sets):
                # Check how much overlap this read has with the current
                # haplotype
                overlap_start = sys.maxsize
                overlap_end = 0
                for la in info.alignments:
                    # a_read is the spanning read by construction (see
                    # `phase` and `get_relevant_read_info`), b_read is the
                    # read in the assembly graph.
                    a_read, b_read = la.get_oriented_reads()
                    if b_read in hap_read_set or b_read in ext_read_set:
                        overlap_start = min(overlap_start, la.arange[0])
                        overlap_end = max(overlap_end, la.arange[1])

                if overlap_start < overlap_end:
                    hap_prob = (overlap_end - overlap_start) / info.overlap_max
                    read_prob += hap_prob

            read_prob /= self.ploidy
            logger.debug("Read probability: %.4f", read_prob)

            if read_prob == 0.0:
                ext_prob += float('-inf')
            else:
                ext_prob += math.log10(read_prob)

        # Calculate P[E] prior, assumed equal for all haplotype sets
        # except for sets with duplicate paths
        multiplicities = Counter(extension)

        multinomial_coeff = 1.0
        multiplicity_sum = 0
        for multiplicity in multiplicities.values():
            multiplicity_sum += multiplicity
            multinomial_coeff *= binom(multiplicity_sum, multiplicity)

        prior = math.log10(multinomial_coeff / num_possible_sets)

        logger.debug("log(P[R|H,E]) = %.3f, log(P[E]) = %.3f", ext_prob, prior)
        logger.debug("log RL[E|H,R] = %.3f", ext_prob + prior)

        return ext_prob + prior
