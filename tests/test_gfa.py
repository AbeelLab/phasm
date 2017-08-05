import tempfile
import logging
import os

import pytest  # noqa

from phasm.assembly_graph import (build_assembly_graph, make_symmetric,
                                  remove_transitive_edges, remove_tips,
                                  clean_graph, remove_diamond_tips,
                                  merge_unambiguous_paths)
from phasm.filter import ContainedReads, MinReadLength, MaxOverhang
from phasm.io import gfa

logging.getLogger().setLevel(logging.INFO)


@pytest.fixture(scope="module")
def assembly_graph():
    path = os.path.join(os.path.dirname(__file__), 'data', 'alignments.gfa')
    reads = {}
    with open(path) as f:
        reads = gfa.gfa2_parse_segments(f)

    filters = [
        ContainedReads(),
        MinReadLength(5000),
        MaxOverhang(1000, 0.8)
    ]

    with open(path) as gfa_file:
        la_iter = map(gfa.gfa2_line_to_la(reads),
                      (l for l in gfa_file if l.startswith('E')))
        la_iter = filter(lambda x: all(f(x) for f in filters), la_iter)

        g = build_assembly_graph(la_iter)

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

        # Free up memory
        del filters

    edges_to_remove = remove_transitive_edges(g, 1000)
    g.remove_edges_from(edges_to_remove)
    make_symmetric(g)

    remove_tips(g, 4, 5000)
    make_symmetric(g)
    clean_graph(g)
    remove_diamond_tips(g)

    remove_tips(g, 4, 5000)
    make_symmetric(g)
    clean_graph(g)
    merge_unambiguous_paths(g)

    return g


@pytest.fixture(scope="module")
def orig_reads():
    path = os.path.join(os.path.dirname(__file__), 'data', 'alignments.gfa')
    with open(path) as f:
        reads = gfa.gfa2_parse_segments(f)
        return reads


def test_gfa_reconstruction(tmpdir, assembly_graph, orig_reads):
    tmpfile = tempfile.NamedTemporaryFile(delete=False)

    with open(tmpfile.name, "w") as f:
        gfa.gfa2_write_graph(f, assembly_graph)

    # Reconstruct
    with open(tmpfile.name) as f:
        graph_reads = gfa.gfa2_parse_segments_with_fragments(f)

    with open(tmpfile.name) as f:
        g = gfa.gfa2_reconstruct_assembly_graph(f, graph_reads, orig_reads)

    for n in assembly_graph.nodes_iter():
        assert n in g

    for u, v, d in assembly_graph.edges_iter(data=True):
        assert u in g
        assert v in g
        assert v in g[u]
        assert d == g[u][v]

    tmpfile.close()
