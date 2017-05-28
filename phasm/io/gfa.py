"""
Utilities to output Graphical Assembly Format-files
===================================================

Tools to convert a NetworkX assembly graph to a GFA file.
"""

from typing import TextIO, Mapping

from phasm.alignments import Read, LocalAlignment, Strand
from phasm.assembly_graph import AssemblyGraph


def gfa_line_to_read(line: str):
    if not line.startswith('S'):
        raise ValueError("Given GFA2 line is not a segment.")

    parts = line.strip().split('\t')
    sid = parts[1].strip()
    length = int(parts[2])

    sequence = None
    if parts[3] != "*":
        sequence = parts[3].encode('ascii')

    return Read(sid, length, sequence)


def _gfa_pos_to_int(pos: str):
    if pos.endswith('$'):
        return int(pos[:-1])

    return int(pos)


def gfa_line_to_la(reads: Mapping[str, Read]):
    def mapper(line: str):
        if not line.startswith('E'):
            raise ValueError('Given GFA2 line is not an edge.')

        parts = line.strip().split('\t')
        sid1 = parts[2].strip()
        sid2 = parts[3].strip()
        arange = tuple(map(_gfa_pos_to_int, parts[4:6]))
        brange = tuple(map(_gfa_pos_to_int, parts[6:8]))
        alignment = parts[8].strip()
        tags = []
        if len(parts) > 9:
            tags.extend(parts[9:])

        a_read = reads[sid1[:-1]]
        b_read = reads[sid2[:-1]]

        if sid1[-1] == '-':
            raise ValueError("A-read as reverse complement, unsupported by "
                             "phasm.")

        strand = Strand.SAME if sid2[-1] == '+' else Strand.OPPOSITE

        return LocalAlignment(a_read, b_read, strand, arange, brange,
                              alignment)

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
