"""
Haplotype phasing using an assembly graph
=========================================

This module contains the algorithm to phase a haplograph to a set of linear DNA
strings, each representing a chromosome copy.
"""

import math
import random
import logging
from itertools import combinations_with_replacement, product
from typing import List, Iterable, Set, Tuple
from collections import OrderedDict, deque, defaultdict, Counter

import networkx
from scipy.stats import norm
from scipy.special import binom

from phasm.alignments import OrientedRead, OrientedDNASegment, MergedReads
from phasm.assembly_graph import AssemblyGraph
from phasm.bubbles import find_superbubbles, superbubble_nodes
from phasm.typing import Node, AlignmentsT, PruneParam, RelevantLA

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

        # Keep the boundaries of each haploblock
        self.block_bounds = []

        if isinstance(copy_from, HaplotypeSet):
            for i in range(ploidy):
                self.haplotypes.append(deque(copy_from.haplotypes[i]))
                self.read_sets.append(set(copy_from.read_sets[i]))
                self.block_bounds.append(list(copy_from.block_bounds[i]))
        else:
            for i in range(ploidy):
                self.haplotypes.append(deque())
                self.read_sets.append(set())
                self.block_bounds.append(list())

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

    def new_block(self):
        for i in range(self.ploidy):
            self.block_bounds.append(len(self.haplotypes[i]))


class BubbleChainPhaser:
    def __init__(self, g: AssemblyGraph,
                 alignments: AlignmentsT,
                 ploidy: int,
                 coverage_model: CoverageModel,
                 min_spanning_reads: int,
                 threshold: PruneParam,
                 prune_factor: PruneParam):
        self.g = g
        self.alignments = alignments
        self.ploidy = ploidy
        self.coverage_model = coverage_model

        self.num_bubbles = 0
        self.prev_bubble = None

        # Keep a record of all reads seen
        self.prev_aligning_reads = set()  # type: Set[OrientedRead]

        # Keep a record of reads playing a role in the assembly graph
        self.prev_graph_reads = set()  # type: Set[OrientedRead]

        self.min_spanning_reads = min_spanning_reads
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
                new_block = self.branch(curr_node, bubbles[curr_node],
                                        bubble_num)
                if not new_block:
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

        # Determine reads that span between paths in previous bubbles and the
        # current bubble.

        # The difference between `cur_bubble_graph_reads` and
        # `cur_bubble_aligning_reads` is that the former one only contains
        # reads that play a role in the assembly graph, and the latter one also
        # contains reads that don't necessarily play a role in the assembly
        # graph, but do align to the current bubble.
        #
        # Additionaly, both variables only include reads that are related
        # to the *interior* of the current bubble, so an alignment to the
        # entrance or exit is ignored. Alignments to the entrance or exit do
        # not provide any information on connecting paths between two
        # different bubbles.
        interior_nodes = (superbubble_nodes(
            self.g, entrance, exit) - {entrance, exit})
        cur_bubble_graph_reads = set(self.get_all_reads(interior_nodes))
        cur_bubble_aligning_reads = self.get_aligning_reads(interior_nodes)
        logger.debug("Current bubble has %d reads used for its interior and "
                     "%d reads aligning to this bubble.",
                     len(cur_bubble_graph_reads),
                     len(cur_bubble_aligning_reads))
        spanning_reads = self.prev_aligning_reads & cur_bubble_aligning_reads

        if (len(spanning_reads) < self.min_spanning_reads and
                len(self.possible_haplotypes) > 1):
            # No reads spanning the previous bubbles and the current one.
            # This means we're actually starting a new haplotype block. But,
            # our `possible_haplotypes` may still be filled with possible
            # candidates, from which none we can link to the current bubble.
            # So, we only keep the most likely candidates, and pick a random
            # one from that list.
            self.prune(1.0, bubble_num)
            random_candidate = random.choice(self.possible_haplotypes)
            random_candidate.new_block()
            self.possible_haplotypes = [random_candidate]

        # Based on the relevant reads, collect all relevant local alignments.
        # We ignore local alignments downstream of the current bubble.
        relevant_la = defaultdict(set)
        total_relevant_la = 0
        for read in spanning_reads:
            for la in self.alignments[read]:
                a_read, b_read = la.get_oriented_reads()

                if (b_read not in cur_bubble_graph_reads and
                        b_read not in self.prev_graph_reads):
                    # Alignment is part of the entrance or exit,
                    # or the alignment does not involve any reads part of the
                    # assembly graph, or the alignment involves a read further
                    # down the graph (beyond the current bubble)
                    continue

                relevant_la[read].add(la)
                total_relevant_la += 1

        if len(spanning_reads) < self.min_spanning_reads:
            logger.info("New haploblock at bubble <%s, %s>, because %s "
                        " spanning reads is not enough.", entrance, exit,
                        len(spanning_reads))
        else:
            logger.info("%d spanning reads inducing %d relevant local "
                        "alignments between bubble <%s, %s> and bubble "
                        "<%s, %s>", len(spanning_reads), total_relevant_la,
                        *self.prev_bubble, entrance, exit)

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
                                        relevant_la, bubble_num)
            )

        self.possible_haplotypes = new_haplotype_sets

        self.prev_bubble = (entrance, exit)
        self.prev_aligning_reads.update(cur_bubble_aligning_reads)
        self.prev_graph_reads.update(cur_bubble_graph_reads)

        return len(spanning_reads) < self.min_spanning_reads

    def generate_new_hsets(self, haplotype_set: HaplotypeSet,
                           possible_paths: List[Tuple[OrientedDNASegment]],
                           relevant_la: RelevantLA,
                           bubble_num: int) -> Iterable[HaplotypeSet]:
        """For a given haplotype set, generate all possible extensions at the
        current bubble. Only yield the new haplotype sets that have a relative
        likelihood above a given threshold."""

        if callable(self.threshold):
            threshold = self.threshold(bubble_num / self.num_bubbles)
        else:
            threshold = self.threshold

        if threshold == 0.0:
            threshold = float('-inf')
        else:
            threshold = math.log10(threshold)

        if bubble_num == 0 or len(relevant_la) < self.min_spanning_reads:
            # For the first bubble the order does not matter, as a permutation
            # in a different order will in the end result in the same haplotype
            # set. The same holds when there are no spanning reads between
            # the previous bubble and the current one. We have no way of
            # linking previous haplotypes to the current one (no read evidence)
            # so this is like starting a new haplotype block
            extension_iter = iter(combinations_with_replacement(possible_paths,
                                                                self.ploidy))
            num_possible_sets = binom(self.ploidy + len(possible_paths) - 1,
                                      self.ploidy)
        else:
            # Otherwise all possible k-tuples of possible paths
            extension_iter = iter(product(possible_paths, repeat=self.ploidy))
            num_possible_sets = self.ploidy**len(possible_paths)

        for extension in extension_iter:
            ext_read_sets = []
            for hap_ext in extension:
                # We index with [1:-1] to ignore the entrance and exit of the
                # bubble
                ext_read_sets.append(set(self.get_all_reads(hap_ext[1:-1])))

            if len(relevant_la) < self.min_spanning_reads:
                new_set = haplotype_set.extend(extension, ext_read_sets)
                yield new_set
            else:
                rl = self.calculate_rl(haplotype_set, extension, ext_read_sets,
                                       relevant_la, num_possible_sets,
                                       bubble_num)

                if rl >= threshold:
                    new_set = haplotype_set.extend(extension, ext_read_sets)
                    new_set.log_rl = rl
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

        if prune_factor > 0.0:
            prune_factor = math.log10(prune_factor)
        else:
            prune_factor = float('-inf')

        max_likelihood = max(hs.log_rl for hs in self.possible_haplotypes)

        logger.info("Max relative likelihood: %.5f", max_likelihood)
        logger.info("Other relative likelihoods: %s", ", ".join(
            "{:.2f}".format(hs.log_rl) for hs in self.possible_haplotypes))

        if max_likelihood != float('-inf'):
            self.possible_haplotypes = [
                hs for hs in self.possible_haplotypes if
                (hs.log_rl - max_likelihood) >= prune_factor
            ]

        num_after = len(self.possible_haplotypes)

        logger.info("Pruned %d unlikely haplotype sets, %d left.",
                    (num_before - num_after), num_after)

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
                     relevant_la: RelevantLA,
                     num_possible_sets: int,
                     bubble_num: int) -> float:
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

        if len(relevant_la) < self.min_spanning_reads:
            return 0.0

        # Calculate P[R|H_set',H_set]
        hs_prob = 0.0
        for read, alignments in relevant_la.items():
            read_prob = 0.0
            for hap_read_set, ext_read_set in zip(hs.read_sets, ext_read_sets):
                num_explained = 0
                for la in alignments:
                    # The a_read is the read that spans two bubbles (by
                    # construction, see also `branch`), so the b_read should
                    # either be a node in the previous haplotype or in the
                    # current extension
                    a_read, b_read = la.get_oriented_reads()
                    if b_read in hap_read_set or b_read in ext_read_set:
                        num_explained += 1

                hap_prob = num_explained / len(alignments)
                read_prob += hap_prob

            read_prob /= self.ploidy

            if read_prob == 0.0:
                hs_prob += float('-inf')
            else:
                hs_prob += math.log10(read_prob)

        # Calculate P[H_set'|H_set] prior, assumed equal for all haplotype sets
        # except for sets with duplicate paths
        multiplicities = Counter(extension)

        multinomial_coeff = 1.0
        multiplicity_sum = 0
        for multiplicity in multiplicities.values():
            multiplicity_sum += multiplicity
            multinomial_coeff *= binom(multiplicity_sum, multiplicity)

        prior = multinomial_coeff / num_possible_sets

        logger.debug("Multiplicities: %s, multinomial coefficient: %d, "
                     "prior: %.3f",
                     ", ".join(map(str, multiplicities.values())),
                     multinomial_coeff, prior)

        hs_prob += math.log10(prior)
        logger.debug("Relative likelihood: %.5f", hs_prob)

        return hs_prob
