"""
Utilities to output Graphical Assembly Format-files
===================================================

Tools to convert a NetworkX assembly graph to a GFA file.
"""

from typing import TextIO

from phasm.assembly_graph import AssemblyGraph


def write_graph(f: TextIO, g: AssemblyGraph, gfa_version: int=1):
    if gfa_version == 1:
        _write_graph_gfa1(f, g)
    else:
        _write_graph_gfa2(f, g)


def _write_graph_gfa1(f: TextIO, g: AssemblyGraph):
    lines = ["H\tVN:Z:1.0\n"]

    segments = set()
    for n in g.nodes_iter():
        segment_name = str(n)[:-1]
        if segment_name in segments:
            continue

        segments.add(segment_name)

        line = "S\t{}\t*\tLN:i:{}\n".format(segment_name, n.length)
        lines.append(line)

    for u, v, d in g.edges_iter(data=True):
        uid = str(u)[:-1]
        u_strand = str(u)[-1:]
        vid = str(v)[:-1]
        v_strand = str(v)[-1:]

        # Fake CIGAR string just indicating overlap length
        overlap = str(d['overlap_len'])+"M"
        line = "L\t{}\t{}\t{}\t{}\t{}\n".format(
            uid, u_strand, vid, v_strand, overlap)
        lines.append(line)

    f.writelines(lines)


def _write_graph_gfa2(f: TextIO, g: AssemblyGraph):
    pass
