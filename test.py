import logging
import argparse
from collections import defaultdict

import networkx

from phasm.io import gfa
from phasm.assembly_graph import (build_assembly_graph,
                                  remove_transitive_edges, clean_graph,
                                  remove_tips, make_symmetric,
                                  merge_unambiguous_paths)
from phasm.filter import ContainedReads, MaxOverhang, MinReadLength
from phasm.walker import build_haplographs

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Build string graph from pairwise overlaps")
    parser.add_argument(
        'gfa_file',
        help="GFA2 file with all reads (segments) and pairwise alignments."
    )

    args = parser.parse_args()
    logger.info("======== STAGE 1: Build Assembly Graph =========")

    logger.info("Pass [1/2] of GFA2 file to import reads (segments)...")
    reads = {}
    with open(args.gfa_file) as f:
        read_iter = map(gfa.parse_gfa_line,
                        (l.strip() for l in f if l.startswith('S')))
        read_iter = map(gfa.gfa_line_to_read, read_iter)
        reads = {r.id: r for r in read_iter}

        logger.info("Read %d reads from GFA2 file.", len(reads))

    logger.info("Pass [2/2] of GFA2 file to import local alignments "
                "(edges)...")
    # logger.setLevel(logging.DEBUG)
    with open(args.gfa_file) as f:
        la_iter = map(gfa.parse_gfa_line,
                      (l.strip() for l in f if l.startswith('E')))
        la_iter = map(gfa.gfa_line_to_la(reads), la_iter)

        # Apply read/alignment filters
        filters = [
            MinReadLength(5000),
            ContainedReads(),
            MaxOverhang(1000, 0.8)
        ]

        la_iter = filter(lambda x: all(f(x) for f in filters), la_iter)

        g = build_assembly_graph(la_iter)
        logger.setLevel(logging.INFO)

        logger.info("Built initial assembly graph with %d nodes and %d "
                    "edges.",
                    networkx.number_of_nodes(g),
                    networkx.number_of_edges(g))

        for f in filters:
            logger.info("Filter %s removed %d alignments.",
                        f.__class__.__name__, f.filtered)

        logger.info("Contained read removal stage 2... (%d contained reads)",
                    len(filters[1].contained_reads))

        for read in filters[1].contained_reads:
            orig = read.with_orientation('-')
            reverse = read.with_orientation('+')

            if orig in g:
                g.remove_node(orig)

            if reverse in g:
                g.remove_node(reverse)

        # Free up memory
        del filters

    logger.info("After contained removal: graph with %d nodes and %d edges.",
                networkx.number_of_nodes(g),
                networkx.number_of_edges(g))

    with open("assembly_graph.gfa", "w") as f:
        gfa.write_graph(f, g)

    networkx.write_graphml(g, "assembly_graph.graphml")

    logger.info("======== STAGE 2: Graph Cleaning =========")

    num_asymm_edges = 0
    edges_to_remove = remove_transitive_edges(g)
    logger.info("Removing %d transitive edges...", len(edges_to_remove))
    g.remove_edges_from(edges_to_remove)
    num_asymm_edges += make_symmetric(g)

    logger.info("Removing tips...")
    num_in_tips, num_out_tips = remove_tips(g)
    num_asymm_edges += make_symmetric(g)

    logger.info("Removing isolated nodes...")
    num_isolated_nodes = clean_graph(g)

    logger.info("Removed %d tip edges, %d isolated nodes, %d asymmetric "
                "edges.",
                num_in_tips+num_out_tips, num_isolated_nodes, num_asymm_edges)

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

    logger.info("Writing graph to file...")

    with open("assembly_graph_reduced.gfa", "w") as f:
        gfa.write_graph(f, g)

    networkx.write_graphml(g, "assembly_graph_reduced.graphml")

    logger.info("======== STAGE 3: Extract and Phase Haplotypes =========")

    # Start with reading all pairwise mappings, so we can easily retreive which
    # reads we need for likelihood calculation
    logger.info("Reading all pairwise mappings...")
    pairwise_mappings = defaultdict(list)
    with open(args.gfa_file) as f:
        la_iter = map(gfa.parse_gfa_line,
                      (l.strip() for l in f if l.startswith('E')))
        la_iter = map(gfa.gfa_line_to_la(reads), la_iter)

        for la in la_iter:
            pairwise_mappings[la.a.id].append(la)

    logger.info("Done.")
    logger.info("Building haplographs...")

    num_components = 0
    num_haplographs = 0
    for i, component in enumerate(
            networkx.weakly_connected_component_subgraphs(g, copy=False)):
        logger.debug("Connected Component %d with %d nodes and %d edges",
                     i, networkx.number_of_nodes(component),
                     networkx.number_of_edges(component))

        for j, haplograph in enumerate(build_haplographs(component)):
            logger.debug("- Haplograph %d with %d nodes and %d edges",
                         j, networkx.number_of_nodes(haplograph),
                         networkx.number_of_edges(haplograph))

            filename = "component{}.haplograph{}.gfa".format(i, j)
            with open(filename, "w") as f:
                gfa.write_graph(f, haplograph)

            num_haplographs += 1

        num_components += 1

    logger.info("Built %d haplographs from %d weakly connected components.",
                num_haplographs, num_components)
