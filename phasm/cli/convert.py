"""
Phasm CLI entry points
"""

import os
import sys
import json
import logging
import argparse

import dinopy

from phasm.io import daligner, gfa
from phasm.io.daligner import Strand
from phasm.utils import random_string

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def fasta2dazzdb(args: argparse.Namespace):
    """Fix the FASTA/FASTQ header/id's to a DAZZ_DB compatible format such that
    these reads can be imported."""

    file_format = args.format
    if not file_format:
        if args.input != sys.stdin:
            filename = args.input.name
            file_ext = filename[filename.rfind('.')+1:]

            file_format = 'fastq' if file_ext in ('fq', 'fastq') else 'fasta'

    if not file_format:
        logger.error("Could not determine file format. Please specify using "
                     "the -f option.")
        return

    if file_format == 'fastq':
        seq_iter = iter(dinopy.FastqReader(args.input).reads(
            quality_values=False))
    else:
        seq_iter = iter(dinopy.FastaReader(args.input).reads(read_names=True))

    if args.input == sys.stdin:
        name = args.name if args.name else random_string(10)
    else:
        name = os.path.basename(args.input.name)

    moviename = daligner.generate_moviename_hash(name)
    name_mapping = {}
    seq_iter = iter(daligner.fix_header(seq_iter, moviename, name_mapping))

    logger.info("Converting FASTA/FASTQ entries...")
    with dinopy.FastaWriter(args.output, force_overwrite=True) as fw:
        fw.write_entries(seq_iter)

    if args.translations:
        logger.info("Writing name mappings to file...")
        json.dump(name_mapping, args.translations)

    logger.info("Done.")


def daligner2gfa(args: argparse.Namespace):
    internal_id_map = {}
    read_lengths = {}
    read_id_trans = {}
    if args.translations:
        if os.path.isfile(args.translations):
            with open(args.translations) as f:
                read_id_trans = json.load(f)
        else:
            logger.warning("Translations file '{}' does not exists, "
                           "ignoring.".format(args.translations_file))

    logger.info("Importing DAZZ_DB reads...")

    args.out.write(gfa.gfa_header(trace_spacing=args.with_trace_points))

    for read in daligner.parse_reads(args.db_input):
        read_id = read['read_id']
        if read_id_trans:
            full_id = daligner.full_id(read)
            read_id = read_id_trans[full_id]

        # If the ID contains white space, only keep the first part
        read_id = read_id.split()[0]

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

        if args.with_trace_points and 'trace_points' in la:
            parts.append(",".join(str(t[1]) for t in la['trace_points']))
        else:
            parts.append("*")

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
        '-T', '--translations', type=str, required=False,
        default=None, metavar='TRANSLATION_FILE',
        help="Give a map which translates DALIGNER sequence ID's to something "
             "else. Should be a JSON file with a single object containing "
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

    fasta2dazzdb_parser = subparsers.add_parser(
        'fasta2dazzdb',
        help="Fix the header (names) of FASTA/FASTQ entries in the given file "
             "to a PacBio compatible format so they can be imported in DAZZ_DB"
             "."
    )
    fasta2dazzdb_parser.set_defaults(func=fasta2dazzdb)

    fasta2dazzdb_parser.add_argument(
        'input', type=argparse.FileType('rb'), default=sys.stdin, nargs='?',
        help="Filename of the FASTA/FASTQ file to read, default stdin"
    )

    fasta2dazzdb_parser.add_argument(
        '-o', '--output', type=argparse.FileType('wb'), default=sys.stdout,
        help="File to write the converted FASTA file to, default stdout"
    )
    fasta2dazzdb_parser.add_argument(
        '-f', '--format', default=None,
        help="Specify the input file format. Can be either fasta or fastq. "
             "Using this setting overrides our own format detection based "
             "on file extension."
    )

    fasta2dazzdb_parser.add_argument(
        '-n', '--name', required=False, default="",
        help="When fixing FASTA headers this script generates a 'moviename' "
             "based on the hash of the filename. When reading from stdin this "
             "filename is not available, and with this option you can specify "
             "what name to use for moviename generation (does not have to be "
             "an actual file). By default it just generates a random string."
    )

    fasta2dazzdb_parser.add_argument(
        '-T', '--translations', default=None, type=argparse.FileType('w'),
        required=False, metavar="TRANSLATIONS_FILE",
        help="Generate a JSON file which contains an new_name: old_name "
             "mapping so you can easility retreive the original name again "
             "for each read."
    )

    args = parser.parse_args()

    if args.func:
        args.func(args)
    else:
        print("No subcommand given, use the", sys.argv[0], "-h to view the "
              "available subcommands")
