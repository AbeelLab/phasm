"""
Phasm CLI entry points
"""

import sys
import json
import logging
import argparse

from phasm.io import daligner, gfa
from phasm.alignments import Strand

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def daligner2gfa(args: argparse.Namespace):
    internal_id_map = {}
    read_lengths = {}
    read_id_trans = {}
    if args.translations:
        read_id_trans = json.load(args.translations)

    logger.info("Importing DAZZ_DB reads...")

    args.out.write(gfa.gfa_header(trace_spacing=args.with_trace_points))

    for read in daligner.parse_reads(args.db_input):
        read_id = read['read_id']
        if read_id_trans:
            full_id = daligner.full_id(read)
            read_id = read_id_trans[full_id]

        # Keep actual read ID to properly define edges below
        internal_id_map[read['read_id']] = read_id

        # Keep read lengths for later
        read_lengths[read['read_id']] = read['length']

        parts = [
            "S", read_id, str(read['length']),
            read['sequence'] if args.with_sequences and 'sequence' in read
            else "*"
        ]

        args.out.write(gfa.gfa_line(*parts))

    logger.info("Reads have been written.")
    logger.info("Importing local alignments...")

    for la in daligner.parse_local_alignments(args.las_input):
        parts = ["E", "*"]

        a_id = internal_id_map[la['a']]
        b_id = internal_id_map[la['b']]

        parts.append(a_id + "+")
        parts.append(b_id + ("+" if la['strand'] == Strand.SAME else "-"))

        arange = list(map(str, la['arange']))
        brange = list(map(str, la['brange']))
        if la['arange'][1] == read_lengths[la['a']]:
            arange[1] += "$"

        if la['brange'][1] == read_lengths[la['b']]:
            brange[1] += "$"

        parts.extend(arange + brange)

        if args.with_trace_points:
            if 'trace_points' in la:
                parts.append(",".join(
                    str(t[1]) for t in la['trace_points']))

        args.out.write(gfa.gfa_line(*parts))

    logger.info("Done. GFA2 file ready.")


def main():
    parser = argparse.ArgumentParser(
        description="A tool to convert several data formats and sources to "
                    "a PHASM compatible input files."
    )

    parser.set_defaults(func=None)
    subparsers = parser.add_subparsers()

    da2gfa_parser = subparsers.add_parser(
        'daligner2gfa', help="Converts your DAZZ_DB and DALIGNER local "
                             "alignments to a GFA2 file."
    )
    da2gfa_parser.set_defaults(func=daligner2gfa)

    da2gfa_parser.add_argument(
        '-s', '--with-sequences', action="store_true", default=False,
        required=False,
        help="Include the full sequence of each read in the resulting GFA file"
    )
    da2gfa_parser.add_argument(
        '-t', '--with-trace-points', type=int, default=None, required=False,
        metavar='SPACING',
        help="To include DALIGNER trace points in the GFA file, set this flag "
             "to the trace point spacing value."
    )
    da2gfa_parser.add_argument(
        '-T', '--translations', type=argparse.FileType('r'), required=False,
        default=None, metavar='TRANSLATION_FILE',
        help="Give a map which translates DALIGNER sequence ID's to something "
             "else. Should be a JSON value with a single object containing "
             "key: value pairs."
    )
    da2gfa_parser.add_argument(
        '-o', '--out', type=argparse.FileType('w'), default=sys.stdout,
        help="Output file, defaults to stdout"
    )

    da2gfa_parser.add_argument(
        'db_input', type=argparse.FileType('r'),
        help="File to read DAZZ_DB DBdump data from. Warning, this should not "
             "be the actual .db but a file containing the DBdump output. See "
             "documentation for recommended usage."
    )

    da2gfa_parser.add_argument(
        'las_input', type=argparse.FileType('r'), default=sys.stdin, nargs='?',
        help="File to read local alignments from (should be in LAdump format)."
             " Defaults to stdin. See documentation for recommended usage."
    )

    args = parser.parse_args()

    if args.func:
        args.func(args)
    else:
        print("No subcommand given, use the", sys.argv[0], "-h to view the "
              "available subcommands")
