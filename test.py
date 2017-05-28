import logging
import argparse
from collections import defaultdict

import networkx
import dinopy

from phasm.io import gfa, sequences
from phasm.assembly_graph import (build_assembly_graph,
                                  remove_transitive_edges, clean_graph,
                                  remove_tips, make_symmetric,
                                  merge_unambiguous_paths)
from phasm.filter import ContainedReads, MaxOverhang, MinReadLength
from phasm.walker import build_haplographs
from phasm.phasing import HaploGraphPhaser, MismatchErrorModel

logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def phase_haplograph(hg, alignments, error_model):
    phaser = HaploGraphPhaser(hg, alignments, 3, error_model, 0.1, 0.8)

    haplotypes = phaser.phase()
    logger.info("Got %d equally likely possible haplotypes", len(haplotypes))

    return haplotypes[0]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Build string graph from pairwise overlaps")
    parser.add_argument(
        'gfa_file',
        help="GFA2 file with all reads (segments) and pairwise alignments."
    )
    parser.add_argument(
        'sequences_file',
        help="Original FASTA/FASTQ file with the sequence data for each read."
    )

    args = parser.parse_args()
    logger.info("======== STAGE 1: Build Assembly Graph =========")

    logger.info("Pass [1/2] of GFA2 file to import reads (segments)...")
    reads = {}
    with open(args.gfa_file) as f:
        read_iter = map(gfa.gfa_line_to_read,
                        (l for l in f if l.startswith('S')))
        reads = {r.id: r for r in read_iter}

        logger.info("Read %d reads from GFA2 file.", len(reads))

    logger.info("Pass [2/2] of GFA2 file to import local alignments "
                "(edges)...")
    with open(args.gfa_file) as f:
        la_iter = map(gfa.gfa_line_to_la(reads),
                      (l for l in f if l.startswith('E')))

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
            filtered = f.filtered

            if f.nodes_to_remove:
                for read in f.nodes_to_remove:
                    orig = read.with_orientation('-')
                    reverse = read.with_orientation('+')

                    if orig in g:
                        filtered += g.degree(orig)
                        g.remove_node(orig)

                    if reverse in g:
                        filtered += g.degree(reverse)
                        g.remove_node(reverse)

            logger.info("Filter %s removed %d alignments.",
                        f.__class__.__name__, f.filtered)

        # Free up memory
        del filters

    logger.info("Final graph: %d nodes and %d edges.",
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
        la_iter = map(gfa.gfa_line_to_la(reads),
                      (l for l in f if l.startswith('E')))

        for la in la_iter:
            a_read, b_read = la.get_oriented_reads()
            pairwise_mappings[a_read].append(la)

    logger.info("Done.")

    sequence_source = sequences.FastaSource(args.sequences_file)
    g.sequence_src = sequence_source

    error_model = MismatchErrorModel(0.01)

    logger.info("Building haplographs...")
    num_components = 0
    num_haplographs = 0

    logger.setLevel(logging.DEBUG)
    with dinopy.FastaWriter("output.fasta", force_overwrite=True) as fw:
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

                filename = "component{}.haplograph{}.graphml".format(i, j)
                networkx.write_graphml(haplograph, filename)

                if networkx.number_of_nodes(haplograph) == 1:
                    seq = g.get_sequence(haplograph.nodes()[0])
                    fw.write_entry(
                        (seq, b"haplotig%d.%d" % (num_haplographs, 0)))
                else:
                    haplotypes = phase_haplograph(haplograph,
                                                  pairwise_mappings,
                                                  error_model)
                    for h, haplotype in enumerate(haplotypes):
                        seq = g.sequence_for_path(haplotype)
                        fw.write_entry(
                            (seq, b"haplotig%d.%d" % (num_haplographs, h)))

                num_haplographs += 1

            num_components += 1

    logger.info("Built %d haplographs from %d weakly connected components.",
                num_haplographs, num_components)
