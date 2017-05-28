import enum
import logging
from typing import Iterable
from collections import defaultdict, OrderedDict

import networkx

from phasm.alignments import (LocalAlignment, AlignmentType, MergedReads,
                              OrientedDNASegment)

logger = logging.getLogger(__name__)


class NodeState(enum.IntEnum):
    VACANT = 0
    IN_PLAY = 1
    ELIMINATED = 2


class AssemblyGraph(networkx.DiGraph):
    adjlist_dict_factory = OrderedDict

    def __init__(self, data: dict=None, **kwargs):
        super().__init__(data, **kwargs)

        self.sequence_src = None

    def subgraph(self, nbunch):
        g = super().subgraph(nbunch)
        g.sequence_src = self.sequence_src

        return g

    @property
    def edge_len(self):
        return self.graph.get('edge_len', 'weight')

    @property
    def overlap_len(self):
        return self.graph.get('overlap_len', 'overlap_len')

    def sort_adjacency_lists(self, reverse=False, weight='weight'):
        for v in self:
            sorted_iter = sorted(self.adj[v].items(),
                                 key=lambda e: e[1][weight])
            if reverse:
                sorted_iter = reversed(sorted_iter)

            self.adj[v] = OrderedDict(sorted_iter)

    def get_sequence(self, read: OrientedDNASegment):
        """Get the sequence for an oriented read. Mostly a convenience function
        which passes the request to `SequenceSource`."""

        if not self.sequence_src:
            raise ValueError("No valid sequence source provided. Cannot get "
                             "DNA sequence of a read.")

        return self.sequence_src.get_sequence(read)

    def sequence_for_path(self, path, edge_len='weight', include_last=True):
        """Get the actual DNA sequence for a path through the assembly graph.
        Requires `self.sequence_src` to be set."""

        if not self.sequence_src:
            raise ValueError("No valid sequence source provided. Cannot spell "
                             "DNA sequence of a path.")

        last = None
        sequence_parts = []
        for u, v, data in path:
            sequence = self.sequence_src.get_sequence(u)
            sequence_parts.append(sequence[:data[edge_len]])

            last = v

        if include_last and last:
            sequence_parts.append(self.sequence_src.get_sequence(last))

        return b"".join(sequence_parts)

    def node_path_edges(self, nodes, data=None):
        """A generator that yields the edges for a path between the given
        nodes. If an edge does not exists in the graph an exception is raised.

        Example::

        >>> g = AssemblyGraph()
        >>> g.add_edges_from([(1, 2), (2, 3), (3, 4)])
        >>> g.node_path_edges([1, 2, 3])
        >>> list(g.node_path_edges(path))
        [(1, 2), (2, 3)]

        """
        if len(nodes) < 2:
            raise ValueError("Not enough nodes to generate a path from")

        node_iter = iter(nodes)

        # Automatically raises StopIteration at the end of the iterator
        node_from = next(node_iter)
        while True:
            node_to = next(node_iter)
            if not self.has_edge(node_from, node_to):
                raise ValueError("Given nodes do not form a path, edge {} "
                                 "not exist".format((node_from, node_to)))
            if data:
                if data is True:
                    yield (node_from, node_to, self[node_from][node_to])
                else:
                    yield (node_from, node_to, self[node_from][node_to][data])
            else:
                yield (node_from, node_to)

            node_from = node_to


def build_assembly_graph(la_iter: Iterable[LocalAlignment],
                         edge_len: str='weight',
                         overlap_len: str='overlap_len') -> AssemblyGraph:
    g = AssemblyGraph(edge_len=edge_len, overlap_len=overlap_len)

    logger.info('Start building assembly graph...')

    for la in la_iter:
        la_type = la.classify()
        a_node, b_node = la.get_oriented_reads()
        a_rev, b_rev = a_node.reverse(), b_node.reverse()

        if la_type == AlignmentType.OVERLAP_AB:
            g.add_edge(a_node, b_node, {
                edge_len: la.arange[0] - la.brange[0],
                overlap_len: la.get_overlap_length()
            })

            g.add_edge(b_rev, a_rev, {
                edge_len: ((len(la.b) - la.brange[1]) -
                           (len(la.a) - la.arange[1])),
                overlap_len: la.get_overlap_length()
            })

            logger.debug('Added edge (%s, %s) with weight %d',
                         a_node, b_node, g[a_node][b_node][edge_len])
            logger.debug('Added edge (%s, %s) with weight %d',
                         b_rev, a_rev, g[b_rev][a_rev][edge_len])
        elif la_type == AlignmentType.OVERLAP_BA:
            g.add_edge(b_node, a_node, {
                edge_len: la.brange[0] - la.arange[0],
                overlap_len: la.get_overlap_length()
            })

            g.add_edge(a_rev, b_rev, {
                edge_len: ((len(la.a) - la.arange[1]) -
                           (len(la.b) - la.brange[1])),
                overlap_len: la.get_overlap_length()
            })

            logger.debug('Added edge (%s, %s) with weight %d',
                         b_node, a_node, g[b_node][a_node][edge_len])
            logger.debug('Added edge (%s, %s) with weight %d',
                         a_rev, b_rev, g[a_rev][b_rev][edge_len])

    logger.info("Built assembly graph with %d nodes and %d edges.",
                networkx.number_of_nodes(g),
                networkx.number_of_edges(g))

    return g


def remove_transitive_edges(g: AssemblyGraph, length_fuzz: int=1000):
    """This function implements the transitive edge reduction algorithm
    described by Myers (2005) [MYERS2005]_.

    It removes edges that are "redundant" in the sense that if there is an
    edge between node `x` and node `y`, an edge between node `y` and node `z`,
    and a node between node `x` and node `z`, the edge between `x` and `z` can
    be removed.

    As a visualisation::

           ___
         /     \    =>   x---y---z
        x---y---z

    .. [MYERS2005] Myers, E. W. (2005). The fragment assembly string graph.
                   Bioinformatics, 21(SUPPL. 2), 79–85.
                   http://doi.org/10.1093/bioinformatics/bti1114
    """

    edges_to_remove = []
    node_state = defaultdict(lambda: NodeState.VACANT)

    logger.info("Transitive edge removal for graph with %d nodes and %d edges",
                networkx.number_of_nodes(g),
                networkx.number_of_edges(g))

    # Ensure that when we iterate over a node's neighbours, we obtain the
    # "shortest" edge first.
    g.sort_adjacency_lists(weight=g.edge_len)
    logger.debug("Sorted adjacency lists")

    num_nodes = networkx.number_of_nodes(g)
    for i, v in enumerate(g):
        v_neighbours = g[v]

        logger.debug("Processing node %s with %d successors (%d/%d)",
                     v, g.out_degree(v), i, num_nodes)
        if not v_neighbours:
            continue

        # Get the last item in the sorted adjacency list, which is the longest
        # edge
        longest_edge = next(reversed(v_neighbours.values()))[g.edge_len]
        longest_edge += length_fuzz

        for w in v_neighbours:
            node_state[w] = NodeState.IN_PLAY

        for w in v_neighbours:
            if node_state[w] == NodeState.IN_PLAY:
                w_neighbours = g[w]
                for x in w_neighbours:
                    if node_state[x] == NodeState.IN_PLAY:
                        total_edge_len = (g[v][w][g.edge_len] +
                                          g[w][x][g.edge_len])
                        if total_edge_len <= longest_edge:
                            node_state[x] = NodeState.ELIMINATED
                            logger.debug("Node %s marked eliminated", x)

        for w in v_neighbours:
            w_neighbours = g[w]
            first = True
            for x in w_neighbours:
                if node_state[x] == NodeState.IN_PLAY:
                    if first:
                        # The first neighbour is the smallest (due to sorting)
                        node_state[x] = NodeState.ELIMINATED
                        logger.debug("Node %s marked eliminated (smallest)", x)

                    if g[w][x][g.edge_len] < length_fuzz:
                        node_state[x] = NodeState.ELIMINATED
                        logger.debug("Node %s marked eliminated", x)

                first = False

        for w in v_neighbours:
            if node_state[w] == NodeState.ELIMINATED:
                edges_to_remove.append((v, w))
                logger.debug("Marked edge (%s, %s) for removal", v, w)
                node_state[w] = NodeState.VACANT

    return edges_to_remove


def remove_outgoing_tips(g: AssemblyGraph, max_tip_len: int=5):
    """Remove short outgoing tips from the assembly graph.

    This function removes short "tips": paths which start at junction or a node
    without incoming edges, ends in a node without any outgoing edges, and with
    a length shorter than `max_tip_len`.

    Example::

                       -> vt1 -> vt2
                     /
        v0 -> v1 -> v2 -> v3 -> v4 -> v5 -- .. -> vn

    The edges (v2, vt1), (vt1, vt2) will be removed.

    Afterwards, it is recommended to delete isolated nodes: nodes without
    incoming or outgoing edges.
    """

    # Clean short tips in outgoing direction
    tips = [n for n in g if g.out_degree(n) == 0]
    num_tip_edges = 0
    for tip in tips:
        is_tip = True

        if g.in_degree(tip) != 1:
            continue

        path = [tip]
        curr_node = tip
        prev = g.predecessors(tip)[0]
        while g.in_degree(curr_node) == 1 and g.out_degree(curr_node) <= 1:
            path.insert(0, prev)

            curr_node = prev
            prev = (g.predecessors(curr_node)[0] if g.in_degree(curr_node) > 0
                    else None)

            logger.debug("Current path: %s", path)
            if len(path) > max_tip_len+1:
                is_tip = False
                break

        if is_tip:
            # Path is a tip, remove it
            logger.debug("Removing tip: %s", path)
            num_tip_edges += len(path)-1

            g.remove_edges_from(g.node_path_edges(path))

    return num_tip_edges


def remove_incoming_tips(g: AssemblyGraph, max_tip_len: int=5):
    """Remove short incoming tips from the assembly graph.

    This function removes short "tips": paths which start at junction or a node
    without incoming edges, ends in a node without any outgoing edges, and with
    a length shorter than `max_tip_len`.

    Example::

                vt1 -> vt2 ----\
                                V
        v0 -> v1 -> v2 -> v3 -> v4 -> v5 -- .. -> vn

    The edges (vt1, vt2), (vt2, v4) will be removed.

    Afterwards, it is recommended to delete isolated nodes: nodes without
    incoming or outgoing edges.
    """

    # Clean short tips in outgoing direction
    tips = [n for n in g if g.in_degree(n) == 0]
    num_tip_edges = 0
    for tip in tips:
        is_tip = True

        if g.out_degree(tip) != 1:
            continue

        path = [tip]
        curr_node = tip
        neighbour = g.neighbors(tip)[0]
        while g.out_degree(curr_node) == 1 and g.in_degree(curr_node) <= 1:
            path.append(neighbour)

            curr_node = neighbour
            neighbour = (g.neighbors(curr_node)[0] if
                         g.out_degree(curr_node) > 0 else None)

            if len(path) > max_tip_len+1:
                is_tip = False
                break

        if is_tip:
            # Path is a tip, remove it
            logger.debug("Removing tip: %s", path)
            num_tip_edges += len(path)-1

            g.remove_edges_from(g.node_path_edges(path))

    return num_tip_edges


def remove_tips(g: AssemblyGraph, max_tip_len: int=3):
    """Remove both small incoming and outgoing tips.

    .. seealso:: remove_incoming_tips, remove_outgoing_tips
    """

    num_incoming_tips = remove_incoming_tips(g, max_tip_len)
    num_outgoing_tips = remove_outgoing_tips(g, max_tip_len)

    return num_incoming_tips, num_outgoing_tips


def remove_short_overlaps(g: AssemblyGraph, drop_ratio: float,
                          sort: bool=False):
    if sort:
        g.sort_adjacency_lists(weight=g.edge_len)

    junction_nodes = (n for n in g.nodes_iter() if g.out_degree(n) > 1)

    for v in junction_nodes:
        v_neighbours = g[v]

        max_ovl = max(w[g.overlap_len] for w in v_neighbours.values())
        _, shortest_edge_target, shortest_edge_ovl = next(
            g.edges_iter(v, data=g.overlap_len))

        if max_ovl != shortest_edge_ovl:
            continue

        threshold = int(round(max_ovl * drop_ratio))

        # Longest edges first (higher chance of short overlap)
        for w, data in reversed(v_neighbours.items()):
            if w == shortest_edge_target:
                # Don't remove the shortest edge
                break

            if data[g.overlap_len] < threshold:
                # Remove this edge
                yield (v, w)
            else:
                break


def make_symmetric(g: AssemblyGraph):
    """This function makes the graph symmetric by removing edges for which
    their reverse complement edge has disappeared. Or in other words, if (u, v)
    is an edge, and (v*, u*) is not an edge, then (u, v) is removed. Here u*
    means the reverse complement of u.

    .. note:: This function does not work if unambiguous paths have been
              merged.
    """

    edges_to_remove = [e for e in g.edges_iter()
                       if not g.has_edge(e[1].reverse(), e[0].reverse())]
    g.remove_edges_from(edges_to_remove)

    return len(edges_to_remove)


def clean_graph(g: AssemblyGraph):
    """Delete nodes that do not have any edges."""

    # Remove nodes without any edges
    isolated_nodes = [n for n in g if g.degree(n) == 0]
    g.remove_nodes_from(isolated_nodes)

    return len(isolated_nodes)


def merge_unambiguous_paths(g: AssemblyGraph):
    """Merge unambiguous (non-branching) paths to a single node.

    This method does not take the reverse complements into account, because
    paths on one strand may be unambiguous, but this does not have to be on the
    reverse strand (because of other reads aligning there). It is therefore not
    recommended to run `make_symmetric` after applying this operation."""

    start_points = [
        n for n in g.nodes_iter() if
        (
            (g.in_degree(n) == 1 and g.out_degree(g.predecessors(n)[0]) > 1) or
            (g.in_degree(n) == 0 or g.in_degree(n) > 1)
        ) and g.out_degree(n) == 1
    ]

    num_merged_nodes = 0

    for start in start_points:
        if g.out_degree(start) != 1:
            continue

        nodes_to_merge = [start]
        curr_node = start
        neighbour = g.neighbors(start)[0]
        while (neighbour and g.out_degree(curr_node) == 1 and
               g.in_degree(neighbour) == 1):
            nodes_to_merge.append(neighbour)

            curr_node = neighbour

            neighbour = (g.neighbors(curr_node)[0] if
                         g.out_degree(neighbour) > 0 else None)

        if len(nodes_to_merge) == 1:
            continue

        logger.debug("Found unambiguous path: %s",
                     ", ".join("({}, in: {}, out: {})".format(
                            n.id, g.in_degree(n), g.out_degree(n)) for n in
                               nodes_to_merge))

        # Create the new node and copy the required edges
        new_id = "[" + "|".join(str(r) for r in nodes_to_merge) + "]"
        prefix_lengths = [l for u, v, l in
                          g.node_path_edges(nodes_to_merge, g.edge_len)]
        new_unmatched_prefix = sum(prefix_lengths)
        new_length = new_unmatched_prefix + len(nodes_to_merge[-1])

        # Keep strand from first node
        new_node = MergedReads(new_id, new_length,
                               nodes_to_merge[0].orientation, nodes_to_merge,
                               prefix_lengths)
        g.add_node(new_node)

        # Incoming edges
        logger.debug("In degree of first node: %d",
                     g.in_degree(nodes_to_merge[0]))
        for u, v, data in g.in_edges_iter(nodes_to_merge[0], data=True):
            g.add_edge(u, new_node, **data)

        # Outgoing edges
        logger.debug("Out degree of last node: %d",
                     g.out_degree(nodes_to_merge[-1]))
        for u, v, data in g.out_edges_iter(nodes_to_merge[-1], data=True):
            # The length of the unmatched prefix probably has changed by
            # merging nodes, so we need to adjust the `edge_len` attribute of
            # all outgoing edges.
            new_data = dict(**data)
            new_data[g.edge_len] += new_unmatched_prefix
            g.add_edge(new_node, v, new_data)

        # Remove merged nodes
        logger.debug("Removing nodes: %s", [str(n) for n in nodes_to_merge])
        g.remove_nodes_from(nodes_to_merge)

        logger.debug("New node in-degree: %d, out-degree: %d.",
                     g.in_degree(new_node), g.out_degree(new_node))

        num_merged_nodes += len(nodes_to_merge)

    return num_merged_nodes
