import sys
import logging
import argparse
import subprocess

import networkx

from phasm.io import daligner, gfa
from phasm.assembly_graph import (build_assembly_graph,
                                  remove_transitive_edges, clean_graph,
                                  remove_tips, remove_short_overlaps,
                                  make_symmetric, merge_unambiguous_paths)
from phasm.filter import ContainedReads, MaxOverhang, MinReadLength

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
        reads = {r.id: r for r in daligner.parse_reads(proc.stdout)}

        logger.info("Read %d reads from DAZZ_DB database.", len(reads))

    la_iter = iter(daligner.parse_local_alignments(reads, args.alignments))

    # Apply read/alignment filters
    filters = [
        MinReadLength(5000),
        ContainedReads(),
        MaxOverhang(1000, 0.8)
    ]

    la_iter = filter(lambda x: all(f(x) for f in filters), la_iter)

    g = build_assembly_graph(la_iter)

    for f in filters:
        logger.info("Filter %s removed %d alignments.", f.__class__.__name__,
                    f.filtered)

    with open("assembly_graph.gfa", "w") as f:
        gfa.write_graph(f, g)

    num_asymm_edges = 0
    edges_to_remove = remove_transitive_edges(g)
    logger.info("Removing %d transitive edges...", len(edges_to_remove))
    g.remove_edges_from(edges_to_remove)
    num_asymm_edges += make_symmetric(g)

    logger.info("Removing tips...")
    num_in_tips, num_out_tips = remove_tips(g)
    num_asymm_edges += make_symmetric(g)

    logger.info("Removing short overlaps...")

    num_short_overlaps = 0
    for ratio in [0.5, 0.6, 0.7]:
        short_overlap_edges = list(remove_short_overlaps(g, ratio))
        g.remove_edges_from(short_overlap_edges)
        num_short_overlaps += len(short_overlap_edges)

    num_asymm_edges += make_symmetric(g)

    logger.info("Removing isolated nodes...")
    num_isolated_nodes = clean_graph(g)

    logger.info("Removed %d tip edges, %d short overlaps, %d isolated nodes, "
                "%d asymmetric edges.",
                num_in_tips+num_out_tips, num_short_overlaps,
                num_isolated_nodes, num_asymm_edges)

    logger.info("Removing tips (stage 2)...")
    num_in_tips, num_out_tips = remove_tips(g)
    num_asymm_edges = make_symmetric(g)
    num_isolated_nodes = clean_graph(g)
    logger.info("Removed %d tip edges, %d isolated nodes, "
                "%d asymmetric edges.", num_in_tips+num_out_tips,
                num_isolated_nodes, num_asymm_edges)

    logger.info("Merging unambiguous paths...")
    num_nodes_merged = merge_unambiguous_paths(g)
    logger.info("Merged %d nodes.", num_nodes_merged)

    logger.info("Done.")
    logger.info("%d/%d nodes have in-degree 1",
                len([n for n in g if g.in_degree(n) == 1]),
                networkx.number_of_nodes(g))
    logger.info("%d/%d nodes have out-degree 1",
                len([n for n in g if g.out_degree(n) == 1]),
                networkx.number_of_nodes(g))

    with open("assembly_graph_reduced.gfa", "w") as f:
        gfa.write_graph(f, g)

    networkx.write_graphml(g, "assembly_graph_reduced.graphml")
