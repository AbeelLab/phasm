"""
This module contains some helper functions to parse data from DAZZ_DB and
DALIGNER.
"""

import hashlib
from typing import Iterable

import numpy

from phasm.alignments import Strand


def generate_moviename_hash(filename: str) -> str:
    """Generate a large integer to be used in the fake PacBio movie name
    based on the SHA256 hash of the filename.
    """

    name_hash = int.from_bytes(
        hashlib.sha256(filename.encode('utf-8')).digest()[:8],
        byteorder='little'
    )

    return str(name_hash)


def fix_header(seq_iter, moviename, name_mapping=None):
    """A generator that "fixes" FASTA headers to a PacBio compatible
    format, to ensure these sequences can be imported in DAZZ_DB.
    """

    seq_id_tpl = "m000000_000000_00000_c{name}/{id}/0_{len}"
    for i, read in enumerate(seq_iter):
        new_name = seq_id_tpl.format(
            name=moviename, id=i, len=len(read.sequence))
        if name_mapping is not None:
            name_mapping[read.name.decode('utf-8')] = new_name

        yield read.sequence, new_name.encode('ascii')


def full_id(read_data: dict) -> str:
    return "{moviename}/{read_id}/{pulse_start}_{pulse_end}".format(
        **read_data)


def parse_reads(input_stream: Iterable[str]) -> Iterable[dict]:
    """Import read metadata from DAZZ_DB.

    Parses read metadata in a DAZZ_DB from an interable `input_stream`,
    which should yield lines of data correesponding to the `DBdump` format
    of DAZZ_DB.

    This function is a generator, and yields each parsed read as a Python
    `dict`.
    """

    current_read_data = {}
    for line in input_stream:
        parts = line.split()
        if line.startswith('R'):
            if len(parts) != 2:
                raise ValueError(
                    "Unexpected input when reading read ID from DAZZ_DB: not "
                    "enough parts (expected 2). Line: '{}'".format(line)
                )

            if current_read_data:
                yield current_read_data
                current_read_data = {}

            current_read_data['read_id'] = parts[1]
        if line.startswith('H'):
            if len(parts) != 3:
                raise ValueError(
                    "Unexpected input when reading reads from DAZZ_DB: not "
                    "enough parts. Line: '{}'".format(line)
                )

            current_read_data['moviename'] = parts[2]
        elif line.startswith('L'):
            if len(parts) != 4:
                raise ValueError(
                    "Unexpected input when reading read lengths from DAZZ_DB: "
                    "not enough parts (expected 4). Line: '{}'".format(line)
                )

            pulse_start, pulse_end = map(int, parts[2:4])
            length = pulse_end - pulse_start
            well = int(parts[1])

            current_read_data['pulse_start'] = pulse_start
            current_read_data['pulse_end'] = pulse_end
            current_read_data['length'] = length
            current_read_data['well'] = well
        elif line.startswith('S'):
            sequence = numpy.array(parts[1].strip().encode('ascii'), dtype='S')
            current_read_data['sequence'] = sequence

    if current_read_data:
        yield current_read_data


def parse_local_alignments(input_stream: Iterable[str]) -> Iterable[dict]:
    """Parse DALIGNER LAdump local alignments.

    This function reads from an iterable `input_stream` the available local
    alignments encoded in DALIGNER LAdump format. In other words, you could
    pipe the output of LAdump to this function, and it yields each local
    alignment with all avaible information as a nice Python `dict`.
    """

    current_la_data = {}
    num_tracepoints = 0
    current_tracepoint = 0

    for line in input_stream:
        parts = line.split()

        # Indicates overlap between two reads, and on which strand
        if line.startswith('P'):
            if (current_la_data and 'a' in current_la_data and 'b' in
                    current_la_data):
                yield current_la_data

                current_la_data = {}
                num_tracepoints = 0
                current_tracepoint = 0

            current_la_data['a'] = parts[1]
            current_la_data['b'] = parts[2]
            current_la_data['strand'] = (Strand.SAME if parts[3] == 'n' else
                                         Strand.OPPOSITE)

        # Indicates alignment range between two reads
        if line.startswith('C'):
            a_start, a_end, b_start, b_end = map(int, parts[1:])
            current_la_data['arange'] = (a_start, a_end)
            current_la_data['brange'] = (b_start, b_end)

        if line.startswith('T'):
            current_la_data['trace_points'] = []
            num_tracepoints = int(parts[1])
            current_tracepoint = 0

        if line.startswith('  '):
            if current_tracepoint >= num_tracepoints:
                raise ValueError(
                    "Received more tracepoints than expected (expected {} "
                    "tracepoints).".format(num_tracepoints)
                )

            current_la_data['trace_points'].append(tuple(map(int, parts)))

        if line.startswith('D'):
            current_la_data['differences'] = int(parts[1])

    if current_la_data and 'a' in current_la_data and 'b' in current_la_data:
        yield current_la_data
        current_la_data = {}
