import os
import sys
import logging
import argparse
from typing import Iterable
from collections import defaultdict

import dinopy
import networkx

from phasm.io import gfa
from phasm.io.sequences import FastaSource
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

    if not args.output:
        args.output = [sys.stdout]

    for f in args.output:
        if f == sys.stdout:
            gfa.write_graph(f, g, args.gfa_version)
        else:
            ext = f.name[f.name.rfind('.')+1:]

            if ext == 'gfa':
                gfa.write_graph(f, g, args.gfa_version)
            elif ext == 'graphml':
                networkx.write_graphml(g, f, encoding='unicode')
            else:
                logger.error("File extension '%s' not recognised, ignoring "
                             "output file %s.", ext, f.name)


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
                    gfa.gfa1_write_graph(f, bubblechain)
                elif args.format == 'gfa2':
                    gfa.gfa2_write_graph(f, bubblechain)
                else:
                    networkx.write_graphml(bubblechain, f, encoding='unicode')

            num_bubblechains += 1

        filename = "component{}.{}".format(i, ext)
        path = os.path.join(args.output_dir, filename)
        with open(path, "w") as f:
            if args.format == 'gfa1':
                gfa.gfa1_write_graph(f, component)
            elif args.format == 'gfa2':
                gfa.gfa2_write_graph(f, component)
            else:
                networkx.write_graphml(component, f, encoding='unicode')

        num_components += 1

    logger.info("Built %d bubblechains from %d weakly connected components.",
                num_bubblechains, num_components)


def phase(args):
    # Original reads are used for assembly graph reconstruction below
    logger.info("Pass [1/2] of GFA2 file to import all original reads "
                "(segments)...")
    with open(args.alignments_gfa) as f:
        reads = gfa.gfa2_parse_segments(f)
    logger.info("Read %d reads from GFA2 file.", len(reads))

    # Original pairwise local alignments are used for likelihood calculation
    logger.info("Pass [2/2] of GFA2 file to import all pairwise local "
                "alignments...")
    read_alignments = defaultdict(set)
    with open(args.alignments_gfa) as f:
        la_iter = map(gfa.gfa2_line_to_la(reads),
                      (l for l in f if l.startswith('E')))

        for la in la_iter:
            a_read, b_read = la.get_oriented_reads()
            read_alignments[a_read].add(la)

    logger.info("Done.")

    # Setup sequence source
    sequence_src = FastaSource(args.reads_fasta)

    with dinopy.FastaWriter(args.output, force_overwrite=True) as fw:
        for bubblechain_gfa in args.bubblechain_gfa:
            logger.info("Bubblechain %s", bubblechain_gfa)
            logger.info("Readig reads and fragments part of assembly graph...")
            with open(bubblechain_gfa) as f:
                graph_reads = gfa.gfa2_parse_segments_with_fragments(f)

            logger.info("Reconstructing assembly graph...")
            with open(bubblechain_gfa) as f:
                g = gfa.gfa2_reconstruct_assembly_graph(f, graph_reads, reads)

            g.sequence_src = sequence_src

            logger.info("Done.")

            logger.info("Start phasing process, ploidy %d...", args.ploidy)
            phaser = BubbleChainPhaser(g, read_alignments, args.ploidy, None,
                                       args.min_spanning_reads, args.threshold,
                                       args.prune_factor)

            for i, (haploblock, include_last) in enumerate(phaser.phase()):
                # Output the DNA sequence for each haplotype
                logger.info("Haploblock %d, building DNA sequences for each "
                            "haplotype...", i)
                id_base = bubblechain_gfa[:bubblechain_gfa.rfind('.')]
                for j, haplotype in enumerate(haploblock.haplotypes):
                    seq = g.sequence_for_path(
                        g.node_path_edges(haplotype, data=True),
                        include_last=include_last
                    )

                    name = "{}.haploblock{}.{}".format(
                        id_base, i, j).encode('ascii')
                    fw.write_entry((seq, name))

            logger.info("Done with %s", bubblechain_gfa)


def main():
    parser = argparse.ArgumentParser(
        description="PHASM: Haplotype-aware de novo genome assembly.")
    parser.set_defaults(func=None)

    parser.add_argument(
        '-v', '--verbose', action='count', default=0, required=False,
        help="Increase verbosity level, number of levels: 0, 1, 2"
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
        '-g', '--gfa-version', choices=(1, 2), default=2,
        help="Which GFA version to use when writing a graph to a GFA file "
             "(default GFA2)."
    )
    layout_io_group.add_argument(
        '-o', '--output', type=argparse.FileType('w'), default=[],
        metavar="FILE", action="append",
        help="Output file (default stdout). If a filename is given, it checks "
             "the file extension for output type. Supported file extensions "
             "are 'graphml' and 'gfa'. This option can be used multiple times "
             "to write multiple files."
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
        '-f', '--format', default="gfa2", choices=("gfa1", "gfa2", "graphml"),
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
    # Phase command
    # ------------------------------------------------------------------------
    phase_parser = subparsers.add_parser(
        'phase',
        help="Phase a bubblechain and output DNA sequences for each haplotype"
             " in FASTA format."
    )
    phase_parser.set_defaults(func=phase)

    phase_parser.add_argument(
        '-p', '--ploidy', type=int, required=True,
        help="The ploidy level."
    )
    phase_parser.add_argument(
        '-t', '--threshold', type=float, default=0.1, required=False,
        help="The minimum relative likelihood of a candidate haplotype set "
             "to be considered for any following bubbles (default: 0.1)."
    )
    phase_parser.add_argument(
        '-d', '--prune-factor', type=float, default=0.5, required=False,
        help="Any candidate haplotype set with a relative likelihood lower "
             "than the given prune factor times the top scoring candidate "
             "will be pruned (default: 0.5)."
    )
    phase_parser.add_argument(
        '-s', '--min-spanning-reads', type=int, default=3, required=False,
        help="If there re less spanning reads between two bubbles than the "
             "given number then PHASM will start a new haploblock."
    )
    phase_parser.add_argument(
        '-o', '--output', type=argparse.FileType('wb'), default=sys.stdout,
        help="Output file (default: stdout)."
    )
    phase_parser.add_argument(
        'reads_fasta',
        help="The FASTA file with all your reads. A FASTA index file should "
             "be present."
    )
    phase_parser.add_argument(
        'alignments_gfa',
        help="The GFA2 file with all pairwise local alignments (used to create"
             " the initial assembly graph). This is a different file than your"
             " bubblechain GFA2 file."
    )
    phase_parser.add_argument(
        'bubblechain_gfa', nargs='+',
        help="The bubblechain graph GFA2 file(s). If given multiple files, "
             "these files will be processed sequentially, but the DNA "
             "sequences will be written to the same file."
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
