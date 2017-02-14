import enum
import logging
from typing import Iterable
from collections import defaultdict, OrderedDict

import networkx

from phasm.alignments import LocalAlignment, Strand, AlignmentType

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


def reverse_id(node_id):
    if node_id[-1] not in ["+", "-"]:
        raise ValueError("An oriented read ID should end with either '+' or "
                         "'-', got '{}'".format(node_id))

    if node_id[-1] == '+':
        return node_id[:-1] + '-'
    else:
        return node_id[:-1] + '+'


def build_assembly_graph(la_iter: Iterable[LocalAlignment]) -> AssemblyGraph:
    g = AssemblyGraph()

    for la in la_iter:
        la_type = la.classify()

        if (la_type == AlignmentType.A_CONTAINED or
                la_type == AlignmentType.B_CONTAINED):
            logger.debug('Contained overlap: %s', la)
            continue

        aid = la.a.id + "+"
        bid = la.b.id + ("+" if la.strand == Strand.SAME else "-")
        r_aid = reverse_id(aid)
        r_bid = reverse_id(bid)

        logger.debug("LA: %s", la)
        logger.debug("%s -> %s, %s -> %s", aid, bid, r_bid, r_aid)

        if la_type == AlignmentType.OVERLAP_AB:
            g.add_edge(aid, bid, {
                'weight': la.arange[0] - la.brange[0]
            })

            g.add_edge(r_bid, r_aid, {
                'weight': ((len(la.b) - la.brange[1]) -
                           (len(la.a) - la.arange[1]))
            })

            logger.debug('Added edge from %s to %s with weight %d',
                         aid, bid, g[aid][bid]['weight'])
            logger.debug('Added edge from %s to %s with weight %d',
                         r_bid, r_aid,
                         g[r_bid][r_aid]['weight'])
        else:
            g.add_edge(bid, aid, {
                'weight': la.brange[0] - la.arange[0]
            })

            g.add_edge(r_aid, r_bid, {
                'weight': ((len(la.a) - la.arange[1]) -
                           (len(la.b) - la.brange[1]))
            })

            logger.debug('Added edge from %s to %s with weight %d',
                         bid, aid, g[bid][aid]['weight'])
            logger.debug('Added edge from %s to %s with weight %d',
                         r_aid, r_bid, g[r_aid][r_bid]['weight'])

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

        longest_edge = max(w[edge_len] for w in v_neighbours.values())
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

    logger.info("Removing %d edges...", len(edges_to_remove))
    for e in edges_to_remove:
        # Not using `g.remove_edges_from` because that function checks if `u`
        # and `v` are existing nodes, which makes it too slow
        u, v = e[:2]
        del g.succ[u][v]
        del g.pred[v][u]

    logger.info("Removing nodes with no edges...")
    g.remove_nodes_from([n for n in g if g.degree(n) == 0])
    logger.info("Done.")
