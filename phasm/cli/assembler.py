import os
import sys
import logging
import argparse
from typing import Iterable, TextIO
from collections import defaultdict

import dinopy
import networkx

from phasm.io import gfa
from phasm.io.sequences import FastaSource
from phasm.typing import ReadMapping, LocalAlignment, AlignmentsT
from phasm.assembly_graph import (build_assembly_graph, clean_graph,
                                  remove_transitive_edges, remove_tips,
                                  make_symmetric, merge_unambiguous_paths,
                                  average_coverage_path, build_bubblechains,
                                  identify_contigs, remove_diamond_tips)
from phasm.filter import (ContainedReads, MaxOverhang, MinReadLength,
                          MinOverlapLength)
from phasm.phasing import BubbleChainPhaser
from phasm.utils import DebugDataLogger

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

    read_alignments = defaultdict(dict)

    def alignment_recorder(la_iter: Iterable[LocalAlignment]) -> Iterable[
            LocalAlignment]:
        nonlocal read_alignments

        for la in la_iter:
            a_read, b_read = la.get_oriented_reads()
            read_alignments[a_read][b_read] = la
            read_alignments[b_read][a_read] = la.switch()

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
    num_in_tips, num_out_tips = remove_tips(g, args.max_tip_length,
                                            args.max_tip_length_bases)
    num_asymm_edges += make_symmetric(g)

    logger.info("Removing isolated nodes...")
    num_isolated_nodes = clean_graph(g)

    logger.info("Removed %d tip edges, %d isolated nodes, %d asymmetric "
                "edges.",
                num_in_tips+num_out_tips, num_isolated_nodes, num_asymm_edges)

    num_diamond_tips = remove_diamond_tips(g)
    logger.info("Removed %d diamond tips", num_diamond_tips)

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


def _write_graphs(g, output_dir, filename_tpl, formats):
    for file_format in formats:
        if file_format.startswith('gfa'):
            version = int(file_format[-1])
            filename = filename_tpl + "gfa"
            path = os.path.join(output_dir, filename)

            with open(path, "w") as f:
                gfa.write_graph(f, g, version)
        else:
            filename = filename_tpl + "graphml"
            path = os.path.join(output_dir, filename)
            with open(path, "w") as f:
                networkx.write_graphml(g, f, encoding='unicode')


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
    num_contigs = 0

    allowed_formats = {'gfa1', 'gfa2', 'graphml'}
    formats = []
    for file_format in args.format.split(','):
        file_format = file_format.strip()

        if file_format in allowed_formats:
            formats.append(file_format)
        else:
            logger.warning("File format '%s' not recognised, ignoring.",
                           file_format)

    if not formats:
        logger.critical("No valid file formats specified.")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    for i, component in enumerate(
            networkx.weakly_connected_component_subgraphs(g, copy=False)):
        logger.info("Connected component %d with %d nodes and %d edges.",
                    i, networkx.number_of_nodes(component),
                    networkx.number_of_edges(component))

        bubblechain_nodes = set()
        for j, bubblechain in enumerate(build_bubblechains(component)):
            logger.info("Found bubblechain #%d with %d nodes and %d edges",
                        j, networkx.number_of_nodes(bubblechain),
                        networkx.number_of_edges(bubblechain))

            for n in bubblechain.nodes_iter():
                component.node[n]['bubblechain'] = j

            bubblechain_nodes.update(bubblechain.nodes_iter())

            filename_tpl = "component{}.bubblechain{}.".format(i, j)
            _write_graphs(bubblechain, args.output_dir, filename_tpl,
                          formats)

            num_bubblechains += 1

        logger.info("Building contigs not part of a bubble chain...")
        for j, path in enumerate(identify_contigs(
                component, bubblechain_nodes, args.min_length)):
            # Update attribute for visualisation in Cytoscape (or other
            # graph viewer)
            if len(path) > 1:
                for u, v in g.node_path_edges(path):
                    component[u][v]['contig'] = j

            filename_tpl = "component{}.contig{}.".format(i, j)
            _write_graphs(g.subgraph(path), args.output_dir, filename_tpl,
                          formats)

            num_contigs += 1

        filename_tpl = "component{}.".format(i)
        _write_graphs(component, args.output_dir, filename_tpl, formats)

        num_components += 1

    logger.info("Built %d bubblechains and %d contigs from %d weakly "
                "connected components.",
                num_bubblechains, num_contigs, num_components)


def _get_read_alignments(f: TextIO, reads: ReadMapping) -> AlignmentsT:
    logger.info("Pass 2 of alignments GFA2 file to import all pairwise local "
                "alignments...")
    read_alignments = defaultdict(dict)

    la_iter = map(gfa.gfa2_line_to_la(reads),
                  (l for l in f if l.startswith('E')))

    for la in la_iter:
        a_read, b_read = la.get_oriented_reads()
        read_alignments[a_read][b_read] = la
        read_alignments[b_read][a_read] = la.switch()

    logger.info("Done.")

    return read_alignments


def phase(args):
    # Original reads are used for assembly graph reconstruction below
    logger.info("Pass 1 of alignments GFA2 file to import all original reads "
                "(segments)...")
    with open(args.alignments_gfa) as f:
        reads = gfa.gfa2_parse_segments(f)
    logger.info("Read %d reads from GFA2 file.", len(reads))

    # Setup sequence source
    sequence_src = FastaSource(args.reads_fasta)
    read_alignments = None

    debug_data_log = DebugDataLogger(args.debug_data)

    with dinopy.FastaWriter(args.output, force_overwrite=True) as fw:
        for gfa_file in args.subgraphs:
            logger.info("Subgraph %s", gfa_file)
            logger.info("Readig reads and fragments part of assembly graph...")
            with open(gfa_file) as f:
                graph_reads = gfa.gfa2_parse_segments_with_fragments(f)

            logger.info("Reconstructing assembly graph...")
            with open(gfa_file) as f:
                g = gfa.gfa2_reconstruct_assembly_graph(f, graph_reads, reads)

            g.sequence_src = sequence_src

            logger.info("Done.")

            logger.info("Start phasing process, ploidy %d...", args.ploidy)
            phaser = BubbleChainPhaser(g, args.ploidy, args.min_spanning_reads,
                                       args.max_bubble_size, args.threshold,
                                       args.prune_factor, args.max_candidates,
                                       args.max_prune_rounds,
                                       args.prune_step_size,
                                       debug_data_log=debug_data_log)

            id_base = os.path.basename(gfa_file[:gfa_file.rfind('.')])
            if len(phaser.bubbles) == 0:
                logger.info("No bubbles found, simple contig path with %d "
                            "nodes.", g.number_of_nodes())
                # This is just a simple "contig" path (linear non-braching
                # path)
                if g.number_of_nodes() == 1:
                    seq = g.get_sequence(g.nodes()[0])
                else:
                    seq = g.sequence_for_path(
                        g.node_path_edges(networkx.topological_sort(g),
                                          data=True),
                        edge_len=g.edge_len
                    )

                fw.write_entry((seq, id_base.encode('ascii')))
            else:
                logger.info("Bubble chain with %d bubbles",
                            len(phaser.bubbles))
                if not read_alignments:
                    with open(args.alignments_gfa) as f:
                        read_alignments = _get_read_alignments(f, reads)

                for i, (haploblock, include_last) in enumerate(
                        phaser.phase(read_alignments)):
                    # Output the DNA sequence for each haplotype
                    logger.info("Haploblock %d, building DNA sequences for "
                                "each haplotype...", i)

                    for j, haplotype in enumerate(haploblock.haplotypes):
                        seq = g.sequence_for_path(
                            g.node_path_edges(haplotype, data=True),
                            include_last=include_last
                        )

                        if haploblock.from_large_bubble:
                            name = "{}.largebubble{}".format(id_base, i)
                        else:
                            name = "{}.haploblock{}.{}".format(id_base, i, j)
                        fw.write_entry((seq, name.encode('utf-8')))

                        if haploblock.from_large_bubble:
                            # Only first haplotype is filled
                            break

                    debug_data_log.haploblock(
                        haploblock, "{}.haploblock{}".format(
                            id_base, i)
                    )

            logger.info("Done with %s", gfa_file)


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
        '-T', '--max-tip-length-bases', type=int, default=5000, required=False,
        help="The maximum length (in bases instead of edges) of a tip "
             "(default 5000)."
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
        '-l', '--min-length', type=int, required=False, default=5000,
        help="Some paths in the assembly graph are not part of a bubble chain,"
             " PHASM will not try to phase these paths but outputs them as "
             "'normal' contigs. With this flag you can specify the minimum "
             "length of a contig to be included (default: 5000)."
    )

    chain_parser.add_argument(
        '-f', '--format', default="gfa2",
        help="Comma separated list of output formats. Supported: gfa1, gfa2, "
             "graphml (default: only GFA2). If multiple formats given, each "
             "bubble chain will get a file in each specified format. This "
             "allows you for example to both write GFA2 and GraphML files "
             "at the same time."
    )
    chain_parser.add_argument(
        '-o', '--output-dir',
        help="Output directory. If the directory does not exist, it will "
             "be created. All identified bubble chains will be written to a "
             "separate file in this directory."
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
        '-s', '--min-spanning-reads', type=int, default=3, required=False,
        help="If there re less spanning reads between two bubbles than the "
             "given number then PHASM will start a new haploblock."
    )
    phase_parser.add_argument(
        '-b', '--max-bubble-size', type=int, default=10, required=False,
        help="The maximum number of simple paths through a bubble. If a "
             "bubble contains more paths from its entrance to exit than the "
             "given number, it is considered too big, and a new haploblock "
             "will be created. The bubble itself will be phased on its own "
             "and not in conjunction with other bubbles. Especially for "
             "larger ploidies you may want to lower this number a bit, as the "
             "number of k-tuples is p^k, where p is the number of paths. "
             "Default value is 10."
    )
    phase_parser.add_argument(
        '-t', '--threshold', type=float, default=1e-3, required=False,
        help="The minimum relative likelihood of a candidate haplotype set "
             "to be considered for any following bubbles (default: 1e-3)."
    )
    phase_parser.add_argument(
        '-d', '--prune-factor', type=float, default=0.1, required=False,
        help="Any candidate haplotype set with a relative likelihood lower "
             "than the given prune factor times the top scoring candidate "
             "will be pruned (default: 0.1)."
    )
    phase_parser.add_argument(
        '-c', '--max-candidates', type=int, default=500, required=False,
        help="At each bubble, limit the number of candidate haplotype sets. "
             "If there more candidates than the given number even after an "
             "initial pruning step, we increasingly prune more stringest "
             "another time, until the number of candidates falls below the "
             "given number (default 500). The maximum number of pruning "
             "rounds can be specified with the option '-r'."
    )
    phase_parser.add_argument(
        '-r', '--max-prune-rounds', type=int, default=9, required=False,
        help="Maximum number of pruning rounds if the number of candidate "
             "haplotype sets is to high (default: 9)."
    )
    phase_parser.add_argument(
        '-S', '--prune-step-size', type=float, default=0.1, required=False,
        help="With each pruning round, increase the prune factor with the "
             "given number (default: 0.1)."
    )

    phase_parser.add_argument(
        '-D', '--debug-data', type=argparse.FileType('w'), default=None,
        required=False,
        help="Output another file containing loads of debug data produced "
             "during the phasing process (optional)."
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
        'subgraphs', nargs='+',
        help="The bubblechain/contig graph GFA2 file(s). If given multiple "
             "files, these files will be processed sequentially, but the DNA "
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
