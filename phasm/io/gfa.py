"""
Utilities to output Graphical Assembly Format-files
===================================================

Tools to convert a NetworkX assembly graph to a GFA file, or reconstructing
a NetworkX graph from a GFA file.

A lot of this module has been written under some time pressure, so it's a bit
of a mess, sorry.
"""

from typing import TextIO, Mapping, NamedTuple, List, Union
from collections import defaultdict
from itertools import zip_longest

from phasm.alignments import Read, LocalAlignment, Strand, MergedReads
from phasm.assembly_graph import AssemblyGraph

_MergedFragment = NamedTuple('MergedFragment', [
    ('id', str), ('length', int), ('reads', List[str]),
    ('prefix_lengths', List[int])
])


class MergedFragment(_MergedFragment):
    def __len__(self):
        return self.length


SegmentMapping = Mapping[str, Union[Read, MergedFragment]]


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


def gfa2_parse_edge(line: str):
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

    return (sid1, sid2, arange, brange, alignment, tags)


def gfa2_line_to_la(reads: Mapping[str, Read]):
    def mapper(line: str):
        if not line.startswith('E'):
            raise ValueError('Given GFA2 line is not an edge.')

        sid1, sid2, arange, brange, alignment, tags = gfa2_parse_edge(line)

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
            fragment_reads = []
            prefix_lengths = []

            for fragment_info in sorted(fragments[segment],
                                        key=lambda elem: elem[2][0]):
                _, external_id, segment_range, fragment_range = fragment_info
                fragment_length = fragment_range[1] - fragment_range[0]
                fragment_reads.append(external_id)
                prefix_lengths.append(fragment_length)

            prefix_lengths.pop()
            reads[segment] = MergedFragment(read.id, length, fragment_reads,
                                            prefix_lengths)

    return reads


def gfa2_reconstruct_assembly_graph(gfa_file: TextIO,
                                    segments: SegmentMapping,
                                    with_orig_reads: Mapping[str, Read]=None,
                                    edge_len: str='weight',
                                    overlap_len: str='overlap_len'):
    """Reconstruct an assembly graph from a GFA2 file. If `with_orig_reads` is
    given, it reconstructs the full assemblygraph with original `OrientedRead`
    and `MergedReads` objects as nodes. Otherwise the nodes will be just
    strings."""

    g = AssemblyGraph()

    def get_orig_oriented_read(oriented_sid: str):
        segment = oriented_sid[:-1]
        return with_orig_reads[segment].with_orientation(oriented_sid[-1])

    nodes_map = {}  # type: Mapping[str, Node]
    # First create node objects, and store them in a dict
    for segment, fragment in segments.items():
        if isinstance(fragment, Read):
            node1 = fragment.with_orientation("+")
            node2 = fragment.with_orientation("-")
        else:
            if with_orig_reads:
                reads = list(map(get_orig_oriented_read, fragment.reads))
            else:
                reads = fragment.reads

            node1 = MergedReads(fragment.id, fragment.length, "+", reads,
                                fragment.prefix_lengths)
            node2 = node1.reverse()

        nodes_map[str(node1)] = node1
        nodes_map[str(node2)] = node2

    # Parse edges, if `with_orig_reads` is given, use the mapping created above
    # to obtain the right nodes.
    line_iter = (l for l in gfa_file if l.startswith('E'))
    for line in line_iter:
        sid1, sid2, arange, brange, alignment, tags = gfa2_parse_edge(line)
        length = arange[0] - brange[0]
        overlap = max(arange[1] - arange[0], brange[1] - brange[0])

        n1 = nodes_map[sid1]
        n2 = nodes_map[sid2]
        g.add_edge(n1, n2, {
            edge_len: length,
            overlap_len: overlap
        })

    return g


def gfa_line(*args) -> str:
    return "\t".join(map(str, args)) + "\n"


def gfa_header(version: str="2.0", trace_spacing: int=None) -> str:
    parts = ["H", "VN:z:{}".format(version)]

    if trace_spacing:
        parts.append("TS:i:{:d}".format(trace_spacing))

    return gfa_line(*parts)


def write_graph(f: TextIO, g: AssemblyGraph, version=2, **kwargs):
    if version == 1:
        return gfa1_write_graph(f, g, **kwargs)
    else:
        return gfa2_write_graph(f, g, **kwargs)


def gfa1_write_graph(f: TextIO, g: AssemblyGraph,
                     with_orig_segments: SegmentMapping=None):
    f.write(gfa_header("1.0"))

    segments = set()
    for n in g.nodes_iter():
        segment_name = str(n)[:-1]
        if segment_name in segments:
            continue

        segments.add(segment_name)

        if with_orig_segments and segment_name in with_orig_segments:
            n = with_orig_segments[segment_name]

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


def gfa2_write_graph(f: TextIO, g: AssemblyGraph):
    f.write(gfa_header("2.0"))

    segments = set()
    for n in g.nodes_iter():
        node_str = str(n)
        segment_name = node_str[:-1]

        if segment_name in segments:
            continue

        parts = ["S", segment_name, len(n), "*"]
        f.write(gfa_line(*parts))
        segments.add(segment_name)

        if isinstance(n, MergedReads):
            # Also output fragment lines to denote which reads are merged
            cur_pos = 0
            total_prefix_lengths = sum(l for l in n.prefix_lengths)

            for read, prefix_len in zip_longest(n.reads, n.prefix_lengths):
                frag_start = 0
                frag_end = (prefix_len if prefix_len else
                            (len(n)-total_prefix_lengths))

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

        parts = ["E", "*", sid1, sid2, a_start, a_end, b_start, b_end, "*"]
        f.write(gfa_line(*parts))
