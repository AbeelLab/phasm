import enum
import logging
from typing import Iterable
from collections import defaultdict, OrderedDict

import networkx

from phasm.alignments import LocalAlignment, AlignmentType

logger = logging.getLogger(__name__)


class NodeState(enum.IntEnum):
    VACANT = 0
    IN_PLAY = 1
    ELIMINATED = 2


class AssemblyGraph(networkx.DiGraph):
    adjlist_dict_factory = OrderedDict

    def sort_adjacency_lists(self, reverse=False, weight='weight'):
        for v in self:
            sorted_iter = sorted(self.adj[v].items(),
                                 key=lambda e: e[1][weight])
            if reverse:
                sorted_iter = reversed(sorted_iter)

            self.adj[v] = OrderedDict(sorted_iter)


def build_assembly_graph(la_iter: Iterable[LocalAlignment]) -> AssemblyGraph:
    g = AssemblyGraph()

    logger.info('Start building assembly graph...')

    for la in la_iter:
        la_type = la.classify()
        a_node, b_node = la.get_oriented_reads()

        if la_type == AlignmentType.OVERLAP_AB:
            g.add_edge(a_node, b_node, {
                'weight': la.arange[0] - la.brange[0],
                'overlap_len': la.get_overlap_length()
            })

            logger.debug('Added edge (%s, %s) with weight %d',
                         a_node, b_node, g[a_node][b_node]['weight'])
        elif la_type == AlignmentType.OVERLAP_BA:
            g.add_edge(b_node, a_node, {
                'weight': la.brange[0] - la.arange[0],
                'overlap_len': la.get_overlap_length()
            })

            logger.debug('Added edge (%s, %s) with weight %d',
                         b_node, a_node, g[b_node][a_node]['weight'])

    logger.info("Built assembly graph with %d nodes and %d edges.",
                networkx.number_of_nodes(g),
                networkx.number_of_edges(g))

    return g


def remove_transitive_edges(g: AssemblyGraph, edge_len: str='weight',
                            length_fuzz: int=10):
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
    g.sort_adjacency_lists(weight=edge_len)
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
        longest_edge = next(reversed(v_neighbours.values()))[edge_len]
        longest_edge += length_fuzz

        for w in v_neighbours:
            node_state[w] = NodeState.IN_PLAY

        for w in v_neighbours:
            if node_state[w] == NodeState.IN_PLAY:
                w_neighbours = g[w]
                for x in w_neighbours:
                    if node_state[x] == NodeState.IN_PLAY:
                        total_edge_len = (g[v][w][edge_len] +
                                          g[w][x][edge_len])
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
                        first = False

                    if g[w][x][edge_len] < length_fuzz:
                        node_state[x] = NodeState.ELIMINATED
                        logger.debug("Node %s marked eliminated", x)

        for w in v_neighbours:
            if node_state[w] == NodeState.ELIMINATED:
                edges_to_remove.append((v, w))
                logger.debug("Marked edge (%s, %s) for removal", v, w)
                node_state[w] = NodeState.VACANT

    return edges_to_remove


def node_path_edges(nodes):
    """
    A generator that yields the edges for a path between the given nodes.

    Example:x

    >>> path = ['n1', 'n2', 'n3']
    >>> list(node_path_edges(path))
    [('n1', 'n2'), ('n2', 'n3')]

    """
    if len(nodes) < 2:
        raise ValueError("Not enough nodes to generate a path from")

    node_iter = iter(nodes)

    # Automatically raises StopIteration at the end of the iterator
    node_from = next(node_iter)
    while True:
        node_to = next(node_iter)
        yield (node_from, node_to)

        node_from = node_to


def remove_tips(g: AssemblyGraph, max_tip_len: int=4):
    """Remove short tips from the assembly graph.

    This function removes short "tips": paths which start at junction or a node
    without incoming edges, ends in a node without any outgoing edges, and with
    a length shorter than `max_tip_len`.

    Example::

                       -- vt1 -- vt2
                     /
        v0 -- v1 -- v2 -- v3 -- v4 -- v5 -- .. -- vn

    The edges (v2, vt1), (vt1, vt2) will be removed.

    Afterwards, it is recommended to delete isolated nodes: nodes without
    incoming or outgoing edges.
    """

    # Clean short tips
    tips = [n for n in g if g.out_degree(n) == 0]
    num_tip_edges = 0
    for tip in tips:
        path = [tip]
        logger.debug("Current path: %s", path)
        while len(path) <= max_tip_len:
            curr_node = path[0]
            pred = g.predecessors(curr_node)

            logger.debug("Predecessors: %s", pred)

            if len(pred) == 1:
                # Exactly one predecessor, extend the path
                path.insert(0, pred[0])

                logger.debug("Predecessor in: %d, out: %d",
                             g.in_degree(pred[0]), g.out_degree(pred[0]))

                if g.out_degree(pred[0]) != 1:
                    # This predecessor is a junction, which means the end of
                    # a tip
                    break
            else:
                # In case of multiple or no predecessors, quit the loop
                break

        if ((g.in_degree(path[0]) == 0 or g.out_degree(path[0]) > 1) and
                len(path) > 1):
            # Path is a tip, remove it
            logger.debug("Removing tip: %s", path)
            num_tip_edges += len(path)-1

            g.remove_edges_from(node_path_edges(path))

    return num_tip_edges


def remove_short_overlaps(g: AssemblyGraph, drop_ratio: float,
                          edge_len: str='weight',
                          overlap_len: str='overlap_len',
                          sort: bool=False):
    if sort:
        g.sort_adjacency_lists(weight=edge_len)

    junction_nodes = (n for n in g.nodes_iter() if g.out_degree(n) > 1)

    for v in junction_nodes:
        v_neighbours = g[v]

        max_ovl = max(w[overlap_len] for w in v_neighbours.values())
        _, shortest_edge_target, shortest_edge_ovl = next(
            g.edges_iter(v, data=overlap_len))

        if max_ovl != shortest_edge_ovl:
            continue

        threshold = int(round(max_ovl * drop_ratio))

        # Longest edges first (higher chance of short overlap)
        for w, data in reversed(v_neighbours.items()):
            if w == shortest_edge_target:
                # Don't remove the shortest edge
                break

            if data[overlap_len] < threshold:
                # Remove this edge
                yield (v, w)
            else:
                break


def clean_graph(g: AssemblyGraph):
    """This function cleans the graph by removing any isolated nodes."""

    # Remove nodes without any edges
    isolated_nodes = [n for n in g if g.degree(n) == 0]
    g.remove_nodes_from(isolated_nodes)

    return len(isolated_nodes)
