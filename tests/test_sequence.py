import os

import pytest  # noqa

from phasm.alignments import Read, MergedReads
from phasm.assembly_graph import AssemblyGraph
from phasm.io.sequences import FastaSource


def test_sequence_for_path():
    path = os.path.join(os.path.dirname(__file__), 'data', 'test.fasta')
    print(path)
    sequence_src = FastaSource(path)
    g = AssemblyGraph()
    g.sequence_src = sequence_src

    r1 = Read("read1", 7).with_orientation('+')
    r2 = Read("read2", 11).with_orientation('+')
    r3 = Read("read3", 15).with_orientation('+')
    r4 = Read("read4", 14).with_orientation('+')

    m1 = MergedReads("merged1", 19, "+", [r2, r3], [4])

    assert g.get_sequence(m1) == b"CATGGGGGGGGACTGCTAC"

    g.add_edge(r1, m1, {'weight': 3, 'overlap_len': 4})
    g.add_edge(m1, r4, {'weight': 14, 'overlap_len': 5})

    seq = g.sequence_for_path(g.node_path_edges([r1, m1, r4], data=True))

    assert seq == b"ATGCATGGGGGGGGACTGCTACGTGTGTGTG"
