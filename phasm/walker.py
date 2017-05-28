"""
Obtain graph representations for chromosomes from the total assembly graph.

This module contains functions to obtain "fused" haplotype contigs. Or in other
words: it helps in obtaining subgraphs from the assembly graph that represent
a set of chromosome copies. This subgraph may contain bubbles and other
non-linear paths that could represents SNP's or larger structural variations
between different chromosome copies.
"""

import logging
from typing import Iterator

import networkx

from phasm.assembly_graph import AssemblyGraph
from phasm.bubbles import find_superbubbles, superbubble_nodes

logger = logging.getLogger(__name__)


def build_haplographs(g: AssemblyGraph,
                      min_nodes: int=1) -> Iterator[AssemblyGraph]:
    # Build dictionary which maps the bubble source to the bubble sink
    logger.info("Searching for superbubbles in the assembly graph...")
    bubbles = {b[0]: b[1] for b in find_superbubbles(g)}
    bubble_sources = set(bubbles.keys())
    bubble_sinks = set(bubbles.values())
    logger.debug("Found superbubbles: %s", bubbles)
    logger.info("Graph has %d superbubbles", len(bubbles))

    # Obtain start nodes, nodes which have no incoming edges or at a junction
    # which is not part of a superbubble.
    # Priority to point without incoming edges
    start_points = [
        n for n in g.nodes_iter() if g.in_degree(n) == 0]

    start_points.extend(n for n in g.nodes_iter() if (
            g.in_degree(n) > 1 and n not in bubble_sources and
            n not in bubble_sinks
        )
    )

    logger.info("Number of start points : %d", len(start_points))

    visited = set()
    for start in start_points:
        if start in visited:
            continue

        subgraph_nodes = set()

        logger.debug("New start point %s", start)
        curr_node = start
        while curr_node:
            if curr_node in bubble_sources:
                logger.debug("Bubble source %s", curr_node)
                # Start of a superbubble, include all nodes of the superbubble
                bubble_nodes = superbubble_nodes(g, curr_node,
                                                 bubbles[curr_node])

                logger.debug("%s is bubble source, this bubble contains %d "
                             "nodes.", curr_node, len(bubble_nodes))
                visited.update(bubble_nodes)
                subgraph_nodes.update(bubble_nodes)

                curr_node = bubbles[curr_node]
            elif curr_node in visited:
                # A bubble sink can also be a bubble source, so that's why we
                # check for a bubble source first above, and afterwards this
                # check if we already visited a given node. If this is not
                # a bubble source and we've seen this node before, then we're
                # looping and we quit building the contig any further.
                logger.debug("%s already visited, stopping", curr_node)
                break
            elif g.out_degree(curr_node) == 1:
                # Simple path to extend
                neighbour = g.neighbors(curr_node)[0]
                visited.add(neighbour)
                subgraph_nodes.add(neighbour)

                curr_node = neighbour
            else:
                # We're either at a node that has no outgoing edges, or a node
                # that is not the source of a superbubble but has multiple
                # outgoing edges. We're not sure what to do now so we quit
                # here.
                logger.debug("Current node %s is a junction with out-degree "
                             "%d", curr_node, g.out_degree(curr_node))
                break

        if len(subgraph_nodes) >= min_nodes:
            yield networkx.subgraph(g, subgraph_nodes)
