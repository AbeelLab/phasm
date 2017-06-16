"""
Utilities to output Graphical Assembly Format-files
===================================================

Tools to convert a NetworkX assembly graph to a GFA file.
"""

from typing import TextIO, Mapping, NamedTuple, List
from collections import defaultdict
from itertools import zip_longest

from phasm.alignments import Read, LocalAlignment, Strand, MergedReads
from phasm.assembly_graph import AssemblyGraph

MergedFragment = NamedTuple('MergedFragment', [
    ('id', str), ('length', int), ('reads', List[str]),
    ('prefix_lengths', List[int])
])


def gfa2_segment_to_read(line: str) -> Read:
    if not line.startswith('S'):
        raise ValueError("Given GFA2 line is not a segment.")

    parts = line.strip().split('\t')
    sid = parts[1].strip()
    length = int(parts[2])

    sequence = None
    if parts[3] != "*":
        sequence = parts[3].encode('ascii').upper()

    return Read(sid, length, sequence)


def gfa2_parse_fragment(line: str):
    if not line.startswith('F'):
        raise ValueError("Given GFA2 line is not a fragment")

    parts = line.strip().split('\t')

    sid = parts[1].strip()
    external_id = parts[2].strip()
    segment_range = tuple(map(_gfa_pos_to_int, parts[3:5]))
    fragment_range = tuple(map(_gfa_pos_to_int, parts[5:7]))

    # Ignore alignment and tags

    return (sid, external_id, segment_range, fragment_range)


def _gfa_pos_to_int(pos: str):
    if pos.endswith('$'):
        return int(pos[:-1])

    return int(pos)


def gfa2_line_to_la(reads: Mapping[str, Read]):
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


def gfa2_parse_segments(f: TextIO):
    read_iter = map(gfa2_segment_to_read, (l for l in f if l.startswith('S')))
    return {r.id: r for r in read_iter}


def gfa2_parse_segments_with_fragments(f: TextIO):
    """Parse all segments and fragments from a GFA2 file and store them in a
    dict.

    In PHASM fragments are used to denote merged reads."""

    segments = {}
    fragments = defaultdict(list)
    for line in f:
        if not line.startswith('S') and not line.startswith('F'):
            continue

        parts = line.strip().split('\t')

        line_type = parts[0].strip()
        segment_name = parts[1].strip()
        if line_type == 'S':
            segments[segment_name] = line

        if line_type == 'F':
            fragments[segment_name].append(gfa2_parse_fragment(line))

    reads = {}
    for segment, line in segments.items():
        if segment not in fragments:
            reads[segment] = gfa2_segment_to_read(line)
        else:
            read = gfa2_segment_to_read(line)
            length = len(read)
            reads = []
            prefix_lengths = []

            for fragment_info in sorted(fragments[segment],
                                        key=lambda elem: elem[2][0]):
                _, external_id, segment_range, fragment_range = fragment_info
                fragment_length = fragment_range[1] - fragment_range[0]
                reads.append(external_id)
                prefix_lengths.append(fragment_length)

            prefix_lengths.pop()
            reads[segment] = MergedFragment(read.id, length, reads,
                                            prefix_lengths)

    return reads


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
    f.write(gfa_header("2.0"))

    segments = set()
    for n in g.nodes_iter():
        segment_name = str(n)[:-1]

        if segment_name in segments:
            continue

        parts = ["S", segment_name, len(n), "*"]
        f.write(gfa_line(*parts))
        segments.add(segment_name)

        if isinstance(n, MergedReads):
            # Also output fragment lines to denote which reads are merged
            cur_pos = 0
            for read, prefix_len in zip_longest(n.reads, n.prefix_lengths):
                frag_start = 0
                frag_end = prefix_len

                segment_start = cur_pos
                segment_end = (cur_pos + prefix_len) if prefix_len else len(n)

                parts = ["F", segment_name, str(read), segment_start,
                         segment_end, frag_start, frag_end, "*"]

                f.write(gfa_line(*parts))
                if prefix_len:
                    cur_pos += prefix_len

    for u, v, d in g.edges_iter(data=True):
        sid1 = str(u)
        sid2 = str(v)

        a_start = d[g.edge_len]
        a_end = len(u)
        b_start = 0
        b_end = d[g.overlap_len]

        parts = ["E", sid1, sid2, a_start, a_end, b_start, b_end, "*"]
        g.write(gfa_line(*parts))
