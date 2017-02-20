import sys
import logging
import argparse
import subprocess

import networkx

from phasm.io import parse_local_alignments, parse_reads
from phasm.assembly_graph import (build_assembly_graph,
                                  remove_transitive_edges, clean_graph)

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Build string graph from pairwise overlaps")
    parser.add_argument('database',
                        help="DAZZ_DB main database file")
    parser.add_argument(
        'alignments', type=argparse.FileType(), nargs='?', default=sys.stdin,
        help="File to read local alignments in LAdump format from"
    )

    args = parser.parse_args()

    reads = {}
    with subprocess.Popen(["DBdump", "-rh", args.database],
                          stdout=subprocess.PIPE,
                          universal_newlines=True) as proc:
        reads = {r.id: r for r in parse_reads(proc.stdout)}

        logger.info("Read %d reads from DAZZ_DB database.", len(reads))

    la_iter = iter(parse_local_alignments(reads, args.alignments))
    g = build_assembly_graph(la_iter)

    networkx.write_graphml(g, 'assembly_graph.graphml')

    edges_to_remove = remove_transitive_edges(g)
    logger.info("Removing %d edges...", len(edges_to_remove))
    g.remove_edges_from(edges_to_remove)

    logger.setLevel(logging.DEBUG)
    logger.info("Cleaning graph... (remove tips, remove isolated nodes)")
    num_tip_edges, num_isolated_nodes = clean_graph(g)

    logger.info("Removed %d tip edges, %d isolated nodes.",
                num_tip_edges, num_isolated_nodes)

    logger.info("Done.")
    logger.info("%d/%d nodes have in-degree 1",
                len([n for n in g if g.in_degree(n) == 1]),
                networkx.number_of_nodes(g))
    logger.info("%d/%d nodes have out-degree 1",
                len([n for n in g if g.out_degree(n) == 1]),
                networkx.number_of_nodes(g))
    networkx.write_graphml(g, 'assembly_graph_reduced.graphml')
