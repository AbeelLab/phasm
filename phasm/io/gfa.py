"""
Utilities to output Graphical Assembly Format-files
===================================================

Tools to convert a NetworkX assembly graph to a GFA file.
"""

from typing import TextIO, Mapping

import numpy
import gfapy

from phasm.alignments import Read, LocalAlignment, Strand
from phasm.assembly_graph import AssemblyGraph


def parse_gfa_line(line: str):
    return gfapy.Line(line)


def gfa_line_to_read(line: gfapy.Line):
    if not str(line).startswith('S'):
        raise ValueError("Given GFA2 line is not a segment.")

    sequence = None
    if not gfapy.is_placeholder(line.sequence):
        sequence = numpy.array(line.sequence.encode('ascii'), dtype='S')

    return Read(line.sid, line.slen, sequence)


def gfa_line_to_la(reads: Mapping[str, Read]):
    def mapper(line: gfapy.Line):
        if not str(line).startswith('E'):
            raise ValueError('Given GFA2 line is not an edge.')

        a_read = reads[line.sid1.name]
        b_read = reads[line.sid2.name]

        if line.sid1.orient == '-':
            raise ValueError("A-read as reverse complement, unsupported by "
                             "phasm.")

        strand = Strand.SAME if line.sid2.orient == '+' else Strand.OPPOSITE

        alignment = None
        if isinstance(line.alignment, gfapy.CIGAR):
            alignment = str(line.alignment)
        elif isinstance(line.alignment, gfapy.Trace):
            alignment = list(line.alignment)

        return LocalAlignment(a_read, b_read, strand,
                              (int(line.beg1), int(line.end1)),
                              (int(line.beg2), int(line.end2)), alignment)

    return mapper


def gfa_line(*args) -> str:
    return "\t".join(args) + "\n"


def gfa_header(version: str="2.0", trace_spacing: int=None) -> str:
    parts = ["H", "VN:z:{}".format(version)]

    if trace_spacing:
        parts.append("TS:i:{:d}".format(trace_spacing))

    return gfa_line(*parts)


def write_graph(f: TextIO, g: AssemblyGraph, gfa_version: int=1):
    if gfa_version == 1:
        _write_graph_gfa1(f, g)
    else:
        _write_graph_gfa2(f, g)


def _write_graph_gfa1(f: TextIO, g: AssemblyGraph):
    f.write(gfa_header("1.0"))

    segments = set()
    for n in g.nodes_iter():
        segment_name = str(n)[:-1]
        if segment_name in segments:
            continue

        segments.add(segment_name)

        parts = ["S", segment_name, "*", "LN:i:{}".format(n.length)]
        f.write(gfa_line(*parts))

    for u, v, d in g.edges_iter(data=True):
        uid = str(u)[:-1]
        u_strand = str(u)[-1:]
        vid = str(v)[:-1]
        v_strand = str(v)[-1:]

        # Fake CIGAR string just indicating overlap length
        overlap = str(d['overlap_len'])+"M"

        parts = ["L", uid, u_strand, vid, v_strand, overlap]
        f.write(gfa_line(*parts))


def _write_graph_gfa2(f: TextIO, g: AssemblyGraph):
    pass
