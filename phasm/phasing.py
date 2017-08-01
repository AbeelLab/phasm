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
from phasm.utils import DebugDataLogger

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
                 max_bubble_size: int,
                 threshold: PruneParam,
                 prune_factor: PruneParam,
                 max_candidates: int,
                 max_prune_rounds: int,
                 prune_step_size: float,
                 debug_data_log: DebugDataLogger=None):
        self.g = g
        self.ploidy = ploidy

        self.start_of_block = True

        # Keep a record of all reads seen
        self.prev_aligning_reads = set()  # type: Set[OrientedRead]

        # Keep a record of reads playing a role in the assembly graph
        self.prev_graph_reads = set()  # type: Set[OrientedRead]

        self.min_spanning_reads = min_spanning_reads
        self.max_bubble_size = max_bubble_size

        self.threshold = threshold
        self.prune_factor = prune_factor
        self.max_candidates = max_candidates
        self.max_prune_rounds = max_prune_rounds
        self.prune_step_size = prune_step_size

        self.candidate_sets = [HaplotypeSet(self.ploidy)]

        self.debug_data_log = (debug_data_log if debug_data_log else
                               DebugDataLogger())

        # Superbubbles without an interior are filtered, to make it a bit
        # easier to identify subgraphs which are just linear paths instead of
        # an actual bubble chain.
        def non_branching_bubble(b):
            """Filter superbubbles with an empty interior (i.e. there are no
            other nodes between the entrance and exit)."""

            s, t = b
            return not (g.out_degree(s) == 1 and g.in_degree(t) == 1 and
                        g.successors(s)[0] == t)

        # find_superbubbles finds superbubbles in reverse topological order
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

    def phase(self, alignments: AlignmentsT) -> Iterable[HaplotypeSet]:
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

            # We're using a list of tuples because tuples are immutable and can
            # be used as dict keys, and more importantly in a
            # `collections.Counter` to determine the multiplicity of path in a
            # later stage.
            possible_paths = list(
                tuple(p) for p in networkx.all_simple_paths(self.g, entrance,
                                                            exit)
            )

            logger.info("%d possible paths through bubble <%r, %r>",
                        len(possible_paths), entrance, exit)

            # Determine reads that span between paths in previous bubbles
            # and the current bubble.
            # The difference between `cur_bubble_graph_reads` and
            # `cur_bubble_aligning_reads` is that the former one only
            # contains reads that play a role in the assembly graph, and
            # the latter one also contains reads that don't necessarily
            # play a role in the assembly graph, but do align to the
            # current bubble.
            #
            # Additionaly, both variables only include reads that are
            # related to the *interior* of the current bubble, so an
            # alignment to the entrance or exit is ignored. Alignments to
            # the entrance or exit do not provide any information on
            # connecting paths between two different bubbles.
            interior_nodes = (superbubble_nodes(
                self.g, entrance, exit) - {entrance, exit})
            cur_bubble_graph_reads = set(self.get_all_reads(
                interior_nodes))
            cur_bubble_aligning_reads = self.get_aligning_reads(
                interior_nodes)
            logger.info("Current bubble has %d reads used for its "
                        "interior and %d reads aligning to this bubble.",
                        len(cur_bubble_graph_reads),
                        len(cur_bubble_aligning_reads))

            # If bubble is too large, use a more simple algorithm to determine
            # the haplotypes
            if len(possible_paths) > self.max_bubble_size:
                logger.info("This bubble is very large, start a new "
                            "haploblock, and use greedy algorithm to get "
                            "haplotypes")
                if not self.start_of_block:
                    yield self.new_block()

                rel_read_info = self.get_relevant_read_info(
                    cur_bubble_aligning_reads, cur_bubble_graph_reads, exit)

                yield self.phase_large_bubble(possible_paths, rel_read_info)
                logger.info("Built haplotypes for large bubble, move to next")

                # Move to the next bubble and continue
                entrance = exit
                bubble_num += 1
                continue

            # Otherwise, just continue as usual and start collecting
            # information for calculating the likelihood
            spanning_reads = (self.prev_aligning_reads &
                              cur_bubble_aligning_reads)

            if (len(spanning_reads) < self.min_spanning_reads and not
                    self.start_of_block):
                # Not enough reads to link the current bubble to the
                # previous bubbles, yield current haplotype block and
                # start a new one.
                logger.info("Starting new haplotype block because we have "
                            "too few spanning reads (%d < %d)",
                            len(spanning_reads),
                            self.min_spanning_reads)
                yield self.new_block()

                logger.info("New block created, start collecting relevant "
                            "reads for this bubble.")

            # Based on the relevant reads, collect all relevant local
            # alignments. We ignore any information downstream of the
            # current bubble.
            relevant_reads = (spanning_reads if not self.start_of_block else
                              cur_bubble_aligning_reads)
            rel_read_info = self.get_relevant_read_info(
                relevant_reads, cur_bubble_graph_reads, exit)

            total_relevant_la = 0
            for read_alignments, _ in rel_read_info.values():
                total_relevant_la += len(read_alignments)

            logger.info("%d relevant reads inducing %d relevant local "
                        "alignments", len(relevant_reads),
                        total_relevant_la)

            self.debug_data_log.new_bubble(entrance, exit, self.start_of_block,
                                           rel_read_info)

            if total_relevant_la == 0:
                logger.info("No relevant alignments, can not phase. Skipping "
                            "this bubble.")

                entrance = exit
                bubble_num += 1
                continue

            new_haplotype_sets = self.branch(possible_paths, rel_read_info)
            new_haplotype_sets = self.prune(new_haplotype_sets,
                                            self.prune_factor)

            if new_haplotype_sets:
                self.candidate_sets = new_haplotype_sets

                # Move to next bubble
                entrance = exit
                self.prev_graph_reads.update(cur_bubble_graph_reads)
                self.prev_aligning_reads.update(cur_bubble_aligning_reads)
                bubble_num += 1
                self.start_of_block = False
            else:
                # None of the branched haplotype sets are above threshold,
                # yield current haploblock and start a new block at the current
                # bubble. We don't move to the next bubble because we want to
                # try this bubble again, but now as start of a new block.
                logger.info("None of the branched candidates at this bubble "
                            "has a relative likelihood above the given "
                            "threshold. Starting a new block. If this happens "
                            "too often consider lowering your threshold.")
                yield self.new_block()

                logger.info("Try again at the same bubble, but now as start "
                            "of a new haplotype block")

        # Only keep the most probable haplotype set(s) if one's available
        if not self.start_of_block:
            self.candidate_sets = self.prune(self.candidate_sets, 1.0)
            logger.info("Got %d equally likely candidate sets, picking a "
                        "random one.", len(self.candidate_sets))
            yield random.choice(self.candidate_sets)

        logger.info("Done")

    def new_block(self):
        # Only keep the most probable candidate
        logger.info("Starting new block, so prune candidates from previous "
                    "bubbles.")
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
                               cur_bubble_graph_reads: Set[OrientedRead],
                               bubble_exit: Node
                               ) -> RelevantReads:
        """Get all alignments of the given spanning reads that either involve
        one of the reads in the interior of one of the previous bubbles or
        one of the reads in the interior of the current bubble.

        This ensures we only obtain alignments that do not involve an entrance
        or exit, and we ignore any alignments in following bubbles.

        Another thing we need to take into account is that a read can span
        beyond the current bubble. This function also determines the largest
        index that does not exceed the current superbubble.
        """

        rel_read_info = {}
        for read in relevant_reads:
            alignments = set()
            for b_read, la in self.alignments[read].items():
                # Ignore any local alignments further down the graph, or local
                # alignments involving the entrance or an exit of a superbubble
                # (graph reads variables only contain interior reads)
                if (b_read in cur_bubble_graph_reads or b_read in
                        self.prev_graph_reads):
                    alignments.add(la)

            # Find the maximum index of this read that does not extend into
            # the next bubble.
            overlap_max = len(read)
            for pred in self.g.predecessors_iter(bubble_exit):
                overlap_max = min(
                    self._get_max_possible_overlap(read, pred, bubble_exit),
                    overlap_max
                )

            rel_read_info[read] = RelevantReadInfo(alignments, overlap_max)

        return rel_read_info

    def _get_max_possible_overlap(self, read: OrientedRead,
                                  pred: Node, exit: Node) -> int:
        # Get the actual edge length from the second to last read to the last
        # read. We need to check if the second to last node is not a merged
        # node.
        edge_length = self.g[pred][exit][self.g.edge_len]
        if isinstance(pred, MergedReads):
            edge_length -= sum(pred.prefix_lengths)

        if pred in self.alignments and read in self.alignments[pred]:
            la = self.alignments[pred][read]
            return edge_length - la.arange[0]
        else:
            return len(read)

    def branch(self, possible_paths: List[Tuple[Node]],
               relevant_reads: RelevantReads) -> List[HaplotypeSet]:
        """Generate all new possible haplotype extensions for a
        bubble (`entrance`, `exit`)."""

        new_haplotype_sets = []
        for haplotype_set in self.candidate_sets:
            new_haplotype_sets.extend(
                self.generate_new_hsets(haplotype_set, possible_paths,
                                        relevant_reads)
            )

        return new_haplotype_sets

    def generate_new_hsets(self, haplotype_set: HaplotypeSet,
                           possible_paths: List[Tuple[OrientedDNASegment]],
                           relevant_reads: RelevantReads
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
        else:
            # Otherwise all possible k-tuples of possible paths, because now
            # order does matter
            extension_iter = iter(product(possible_paths, repeat=self.ploidy))

        num_possible_sets = len(possible_paths)**self.ploidy

        for extension in extension_iter:
            ext_read_sets = []
            # Get graph reads of the extension
            for hap_ext in extension:
                # We index with [1:-1] to ignore the entrance and exit of the
                # bubble
                ext_read_sets.append(set(self.get_all_reads(hap_ext[1:-1])))

            rl = self.calculate_rl(haplotype_set, extension, ext_read_sets,
                                   relevant_reads, num_possible_sets)

            if rl >= threshold:
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

        log_prune_factor = math.log10(prune_factor)
        max_likelihood = max(hs.log_rl for hs in candidate_sets)

        logger.info("Max log relative likelihood: %.5f", max_likelihood)
        logger.debug("Other log relative likelihoods: %s", ", ".join(
            "{:.2f}".format(hs.log_rl) for hs in candidate_sets))

        if max_likelihood != float('-inf'):
            candidate_sets = [
                hs for hs in candidate_sets if
                (hs.log_rl - max_likelihood) >= log_prune_factor
            ]

            num_after = len(candidate_sets)
            logger.info("Pruned %d unlikely haplotype sets with prune factor "
                        "%g, %d left.", (num_before - num_after), prune_factor,
                        num_after)

            prune_rounds = 1
            while (num_after > self.max_candidates and prune_rounds
                   < self.max_prune_rounds and prune_factor < 1.0):
                prune_factor += self.prune_step_size
                log_prune_factor = math.log10(prune_factor)

                num_before = len(candidate_sets)

                candidate_sets = [
                    hs for hs in candidate_sets if
                    (hs.log_rl - max_likelihood) >= log_prune_factor
                ]

                prune_rounds += 1
                num_after = len(candidate_sets)

                logger.info("Additional pruning round with new prune factor %g"
                            " removed %d candidates, %d candidates left.",
                            prune_factor, (num_before - num_after), num_after)
        else:
            logger.info("Can't prune if max likelihood = -inf")

        return candidate_sets

    def get_aligning_reads(self, nodes: List[Node]) -> Set[OrientedRead]:
        aligning_reads = set()  # type: Set[OrientedRead]

        for read in self.get_all_reads(nodes):
            aligning_reads.add(read)
            if read in self.alignments:
                for b_read, la in self.alignments[read].items():
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

        if len(relevant_reads) < self.min_spanning_reads:
            return float('-inf')

        logger.debug("Haplotype read sets: %s", hs.read_sets)
        logger.debug("Extension read sets: %s", ext_read_sets)

        # Calculate P[SR|H,E]
        p_sr = 0.0
        read_probs = []
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
                        overlap_end = max(overlap_end, min(
                            info.overlap_max, la.arange[1]))

                if overlap_start < overlap_end:
                    hap_prob = (overlap_end - overlap_start) / info.overlap_max
                    read_prob += hap_prob

            read_prob /= self.ploidy
            read_probs.append(read_prob)
            logger.debug("Read probability: %.4f", read_prob)
            p_sr += read_prob

        p_sr /= len(relevant_reads)

        if p_sr == 0.0:
            p_sr = float('-inf')
        else:
            p_sr = math.log10(p_sr)

        # Calculate P[E] prior, assumed equal for all haplotype sets
        # except for sets with duplicate paths
        multiplicities = Counter(extension)

        multinomial_coeff = 1.0
        multiplicity_sum = 0
        for multiplicity in multiplicities.values():
            multiplicity_sum += multiplicity
            multinomial_coeff *= binom(multiplicity_sum, multiplicity)

        prior = math.log10(multinomial_coeff / num_possible_sets)

        logger.debug("log P[R|H,E] = %.3f, log P[E] = %.3f", p_sr, prior)
        logger.debug("log RL[E|H,R] = %.3f", p_sr + prior)

        self.debug_data_log.candidate_set(hs, extension, read_probs, p_sr,
                                          prior)

        return p_sr + prior

    def phase_large_bubble(self, possible_paths: List[Tuple[Node]],
                           relevant_reads: RelevantReads):
        """Large superbubbles have too many possible configurations to check
        them all. Therefore we implement a simple greedy algorithm that just
        outputs the top `k` best supported paths by read data through this
        superbubble."""

        logger.info("Calculating score for each path through this bubble...")
        scored_paths = {}
        path_read_sets = {}
        for i, path in enumerate(possible_paths):
            # Ignore entrance and exit
            path_read_set = set(self.get_all_reads(path[1:-1]))
            path_read_sets[path] = path_read_set

            path_prob = 0.0
            for read, info in relevant_reads.items():
                # Check how much overlap this read has with the current
                # path
                overlap_start = sys.maxsize
                overlap_end = 0
                for la in info.alignments:
                    # a_read is the spanning read by construction (see
                    # `phase` and `get_relevant_read_info`), b_read is the
                    # read in the path
                    a_read, b_read = la.get_oriented_reads()
                    if b_read in path_read_set:
                        overlap_start = min(overlap_start, la.arange[0])
                        overlap_end = max(overlap_end, min(
                            info.overlap_max, la.arange[1]))

                if overlap_start < overlap_end:
                    read_prob = (
                        (overlap_end - overlap_start) / info.overlap_max)
                    path_prob += read_prob

            path_prob /= len(relevant_reads)
            logger.debug("Path %/%d, score: %g", i+1, len(possible_paths),
                         path_prob)
            scored_paths[path] = path_prob

        sorted_paths = list(sorted(scored_paths, key=scored_paths.get,
                                   reverse=True))
        sorted_paths = sorted_paths[:self.ploidy]

        hs = HaplotypeSet(self.ploidy)
        for i, path in enumerate(sorted_paths):
            hs.haplotypes[i].extend(path)
            hs.read_sets[i].update(path_read_sets[path])
        logger.info("Done.")

        return hs
