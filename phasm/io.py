"""
This module contains some helper functions to parse data from other sources.
"""

from typing import Mapping, Iterable

from phasm.alignments import Read, LocalAlignment, Strand


def parse_reads(input_stream: Iterable[str]) -> Iterable[Read]:
    """Import read metadata from DAZZ_DB.

    Parses read metadata in a DAZZ_DB from an interable `input_stream`,
    which should yield lines of data correesponding to the `DBdump` format
    of DAZZ_DB.

    This function is a generator, and yields each parsed read.
    """

    moviename = ""
    read_id = ""
    for line in input_stream:
        parts = line.split()
        if line.startswith('R'):
            if len(parts) != 2:
                raise ValueError(
                    "Unexpected input when reading read ID from DAZZ_DB: not "
                    "enough parts (expected 2). Line: '{}'".format(line)
                )

            read_id = parts[1]
        if line.startswith('H'):
            if len(parts) != 3:
                raise ValueError(
                    "Unexpected input when reading reads from DAZZ_DB: not "
                    "enough parts. Line: '{}'".format(line)
                )

            moviename = parts[2]
        elif line.startswith('L'):
            if len(parts) != 4:
                raise ValueError(
                    "Unexpected input when reading read lengths from DAZZ_DB: "
                    "not enough parts (expected 4). Line: '{}'".format(line)
                )

            if not read_id or not moviename:
                raise ValueError(
                    "Unexpected output from DBdump, got a length (L) line "
                    "before the read (R) and header (H) lines."
                )

            pulse_start, pulse_end = map(int, parts[2:4])
            length = pulse_end - pulse_start

            read = Read(read_id, moviename, int(parts[1]), pulse_start,
                        pulse_end, length)
            yield read


def parse_local_alignments(reads: Mapping[str, Read],
                           input_stream: Iterable[str]) -> \
                           Iterable[LocalAlignment]:
    """Parse DALIGNER LAdump local alignments.

    This function reads from an iterable `input_stream` the available local
    alignments encoded in DALIGNER LAdump format. In other words, you could
    pipe the output of LAdump to this function, and it yields each local
    alignment with all avaible information as a nice Python `namedtuple`
    (`LocalAlignment`).
    """

    a = b = strand = a_range = b_range = differences = trace_points = None
    num_tracepoints = 0
    current_tracepoint = 0

    for line in input_stream:
        parts = line.split()

        # Indicate overlap between two reads, and on which strand
        if line.startswith('P'):
            if a and b:
                alignment = LocalAlignment(
                    a, b, strand, a_range, b_range, differences, trace_points
                )

                yield alignment

                a = b = strand = a_range = b_range = differences = None
                trace_points = None

            a = reads[parts[1]]
            b = reads[parts[2]]
            strand = Strand.SAME if parts[3] == 'n' else Strand.OPPOSITE

        # Indicate alignment range between two reads
        if line.startswith('C'):
            a_start, a_end, b_start, b_end = map(int, parts[1:])
            a_range = (a_start, a_end)
            b_range = (b_start, b_end)

        if line.startswith('T'):
            trace_points = []
            num_tracepoints = int(parts[1])
            current_tracepoint = 0

        if line.startswith('  '):
            if current_tracepoint >= num_tracepoints:
                raise ValueError(
                    "Received more tracepoints than expected (expected {} "
                    "tracepoints).".format(num_tracepoints)
                )

            trace_points.append(tuple(map(int, parts)))

        if line.startswith('D'):
            differences = int(parts[1])
