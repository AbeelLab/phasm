import os
import sys
import logging
import argparse
from typing import Iterable
from collections import defaultdict

import networkx

from phasm.io import gfa
from phasm.typing import LocalAlignment
from phasm.assembly_graph import (build_assembly_graph, clean_graph,
                                  remove_transitive_edges, remove_tips,
                                  make_symmetric, merge_unambiguous_paths,
                                  average_coverage_path, build_bubblechains)
from phasm.filter import (ContainedReads, MaxOverhang, MinReadLength,
                          MinOverlapLength)
from phasm.phasing import BubbleChainPhaser

logger = logging.getLogger()


def layout(args):
    logger.info("======== STAGE 1: Build Assembly Graph =========")

    logger.info("Pass [1/2] of GFA2 file to import reads (segments)...")
    reads = {}
    with open(args.gfa_file) as f:
        reads = gfa.gfa2_parse_segments(f)
    logger.info("Read %d reads from GFA2 file.", len(reads))

    logger.info("Pass [2/2] of GFA2 file to import local alignments "
                "and build assembly graph...")

    read_alignments = defaultdict(set)

    def alignment_recorder(la_iter: Iterable[LocalAlignment]) -> Iterable[
            LocalAlignment]:
        nonlocal read_alignments

        for la in la_iter:
            a_read, b_read = la.get_oriented_reads()
            read_alignments[a_read].add(la)

            yield la

    filters = [ContainedReads()]

    if args.min_read_length:
        filters.append(MinReadLength(args.min_read_length))

    if args.min_overlap_length:
        filters.append(MinOverlapLength(args.min_overlap_length))

    filters.append(MaxOverhang(args.max_overhang_abs,
                               args.max_overhang_rel))

    with open(args.gfa_file) as gfa_file:
        la_iter = map(gfa.gfa2_line_to_la(reads),
                      (l for l in gfa_file if l.startswith('E')))
        la_iter = alignment_recorder(la_iter)
        la_iter = filter(lambda x: all(f(x) for f in filters), la_iter)

        g = build_assembly_graph(la_iter)

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

    logger.info("======== STAGE 2: Graph Cleaning =========")

    num_asymm_edges = 0
    edges_to_remove = remove_transitive_edges(g, args.length_fuzz)
    logger.info("Removing %d transitive edges...", len(edges_to_remove))
    g.remove_edges_from(edges_to_remove)
    num_asymm_edges += make_symmetric(g)

    logger.info("Removing tips...")
    num_in_tips, num_out_tips = remove_tips(g, args.max_tip_length)
    num_asymm_edges += make_symmetric(g)

    logger.info("Removing isolated nodes...")
    num_isolated_nodes = clean_graph(g)

    logger.info("Removed %d tip edges, %d isolated nodes, %d asymmetric "
                "edges.",
                num_in_tips+num_out_tips, num_isolated_nodes, num_asymm_edges)

    logger.info("Removing tips (stage 2)...")
    num_in_tips, num_out_tips = remove_tips(g, args.max_tip_length)
    num_asymm_edges = make_symmetric(g)
    num_isolated_nodes = clean_graph(g)
    logger.info("Removed %d tip edges, %d isolated nodes, "
                "%d asymmetric edges.", num_in_tips+num_out_tips,
                num_isolated_nodes, num_asymm_edges)

    logger.info("Merging unambiguous paths...")
    num_nodes_merged = merge_unambiguous_paths(g)
    logger.info("Merged %d nodes.", num_nodes_merged)

    logger.info("Done.")

    logger.info("Calculating average coverage for each edge...")
    for u, v in g.edges_iter():
        g[u][v]['avg_coverage'] = average_coverage_path(g, read_alignments,
                                                        [u, v])

    logger.info("Writing graph...")

    if args.format == 'gfa1':
        gfa.gfa1_write_graph(args.output, g)
    elif args.format == 'gfa2':
        gfa.gfa2_write_graph(args.output, g)
    elif args.format == 'graphml':
        networkx.write_graphml(g, args.output, encoding='unicode')
    else:
        logger.critical("Invalid output format specified: %s", args.format)
        sys.exit(1)


def chain(args):
    logger.info("Readig reads and fragments part of the assembly graph...")
    graph_reads = {}
    with open(args.graph_gfa) as f:
        graph_reads = gfa.gfa2_parse_segments_with_fragments(f)

    logger.info("Reconstructing assembly graph...")

    with open(args.graph_gfa) as f:
        g = gfa.gfa2_reconstruct_assembly_graph(f, graph_reads)

    logger.info("Done.")
    logger.info("Enumerate weakly connected components in the graph...")
    num_components = 0
    num_bubblechains = 0

    if args.format.startswith("gfa"):
        ext = "gfa"
    elif args.format == "graphml":
        ext = "graphml"
    else:
        logger.critical("Invalid output format specifified: %s", args.format)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    for i, component in enumerate(
            networkx.weakly_connected_component_subgraphs(g, copy=False)):
        logger.info("Connected component %d with %d nodes and %d edges.",
                    i, networkx.number_of_nodes(component),
                    networkx.number_of_edges(component))

        for j, bubblechain in enumerate(build_bubblechains(component)):
            logger.info("Found bubblechain #%d with %d nodes and %d edges",
                        j, networkx.number_of_nodes(bubblechain),
                        networkx.number_of_edges(bubblechain))

            for n in bubblechain.nodes_iter():
                component.node[n]['bubblechain'] = j

            filename = "component{}.bubblechain{}.{}".format(i, j, ext)
            path = os.path.join(args.output_dir, filename)
            with open(path, "w") as f:
                if args.format == 'gfa1':
                    gfa.gfa1_write_graph(f, bubblechain, graph_reads)
                elif args.format == 'gfa2':
                    gfa.gfa2_write_graph(f, bubblechain, graph_reads)
                else:
                    networkx.write_graphml(bubblechain, f, encoding='unicode')

            num_bubblechains += 1

        filename = "component{}.{}".format(i, ext)
        path = os.path.join(args.output_dir, filename)
        with open(path, "w") as f:
            if args.format == 'gfa1':
                gfa.gfa1_write_graph(f, component, graph_reads)
            elif args.format == 'gfa2':
                gfa.gfa2_write_graph(f, component, graph_reads)
            else:
                networkx.write_graphml(component, f, encoding='unicode')

        num_components += 1

    logger.info("Built %d bubblechains from %d weakly connected components.",
                num_bubblechains, num_components)


def main():
    parser = argparse.ArgumentParser(
        description="PHASM: Haplotype-aware de novo genome assembly.")
    parser.set_defaults(func=None)

    parser.add_argument(
        '-v', '--verbose', action='count', default=0, required=False,
        help="Increase verbosity level with each usage."
    )

    subparsers = parser.add_subparsers()

    # -------------------------------------------------------------------------
    # Layout command
    # -------------------------------------------------------------------------

    layout_parser = subparsers.add_parser(
        'layout', help="Build an assembly graph based on pairwise local "
                       "alignments."
    )
    layout_parser.set_defaults(func=layout)

    alignment_group = layout_parser.add_argument_group(
        'Alignment and read filtering')
    graph_cleaning_group = layout_parser.add_argument_group('Graph cleaning')
    layout_io_group = layout_parser.add_argument_group('Input/output')

    alignment_group.add_argument(
        '-l', '--min-read-length', type=int, required=False, default=0,
        metavar="LENGTH",
        help="Filter reads smaller than the given length (default: disabled)"
    )
    alignment_group.add_argument(
        '-s', '--min-overlap-length', type=int, required=False, default=0,
        metavar="LENGTH",
        help="Minimum length of the overlap between two reads, otherwise "
             "this alignment is ignored. Default is disabled, because it's "
             "something that's usually handled by your overlapper."
    )
    alignment_group.add_argument(
        '-a', '--max-overhang-abs', type=int, default=1000, required=False,
        metavar="LENGTH",
        help="Max absolute overhang length (default: 1000)."
    )
    alignment_group.add_argument(
        '-r', '--max-overhang-rel', type=float, default=0.8, required=False,
        metavar="FRACTION",
        help="Max overhang length as fraction of the overlap length (default: "
             "0.8)."
    )
    graph_cleaning_group.add_argument(
        '-t', '--max-tip-length', type=int, default=4, required=False,
        metavar="NUM",
        help="Maximum number of edges of a path to be called a tip "
             "(default: 4)."
    )
    graph_cleaning_group.add_argument(
        '-F', '--length-fuzz', type=int, default=1000, required=False,
        metavar="LENGTH",
        help="Transitive reduction length fuzz parameter (default: 1000). "
             "See Myers (2005). "
    )

    layout_io_group.add_argument(
        '-f', '--format', choices=('gfa1', 'gfa2', 'graphml'), default='gfa2',
        help="Output format (default gfa2)."
    )
    layout_io_group.add_argument(
        '-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
        metavar="FILE",
        help="Output file (default stdout)"
    )
    layout_io_group.add_argument(
        'gfa_file', help="Input GFA2 file with all pairwise local alignments."
    )

    # ------------------------------------------------------------------------
    # Chain command
    # ------------------------------------------------------------------------
    chain_parser = subparsers.add_parser(
        'chain', help="Identify and build bubblechains")
    chain_parser.set_defaults(func=chain)

    chain_parser.add_argument(
        '-f', '--format', default="gfa2",
        help="Output format (default: GFA2)."
    )
    chain_parser.add_argument(
        '-o', '--output-dir',
        help="Output directory. If the directory does not exist, it will "
             "be created."
    )
    chain_parser.add_argument(
        'graph_gfa',
        help="The assembly graph in GFA2 format. Other graph formats are not "
             "supported. Note that this is a different file than the GFA2 file"
             " with pairwise local alignments."
    )

    # ------------------------------------------------------------------------
    # Argument parsing
    # ------------------------------------------------------------------------
    args = parser.parse_args()

    # Setup logging
    logger.setLevel(logging.INFO)
    phasm_logger = logging.getLogger('phasm')
    phasm_logger.setLevel(logging.WARNING)

    formatter = logging.Formatter("{asctime} - {levelname}:{name}:{message}",
                                  style="{")
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if args.verbose > 0:
        phasm_logger.setLevel(logging.INFO)

    if args.verbose > 1:
        logger.setLevel(logging.DEBUG)
        phasm_logger.setLevel(logging.DEBUG)

    if not args.func:
        parser.print_help()
    else:
        args.func(args)
