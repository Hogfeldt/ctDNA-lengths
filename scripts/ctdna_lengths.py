# Name: ctDNA-lengths
#
# Description: A script for extracting the length distribution
#              of ctDNA WGS data. For more information and
#              examples of how to run the script follow the url
#              below.
#
# Author (primary): Per Høgfeldt
#
# Affiliation: Department of Molecular Medicine,
#              Aarhus University Hospital,
#              Denmark
#
# Url: https://github.com/hogfeldt/...
#
# License: GNU General Public License v3 (GPLv3)
#
###########################################
# --------------| Imports |---------------#
###########################################
# Python stdlib
import argparse
import os
import sys
import logging
import re
import math
import base64
from datetime import datetime
from itertools import cycle
from typing import NamedTuple

# Numeric packages
import numpy as np

# Bioinformatic packages
import pysam

###########################################
# ------------| BED Parser |--------------#
###########################################


class Token(NamedTuple):
    type: str
    value: str
    line: int
    column: int


def tokenize(text):
    # Compile bed specification tokens to regexp
    token_specification = [
        ("CHROMOSOME", r"chr[3-9]|chr1[0-9]?|chr2[0-2]?|chrX|chrY"),
        ("NUMBER", r"\d+"),
        ("TAB", r"\t"),
        ("NEWLINE", r"\n"),
        ("COMMENT", r"#.*\n"),
        ("ANY", r"."),
    ]
    tok_regex = "|".join("(?P<%s>%s)" % pair for pair in token_specification)

    # Set variables
    line_num = 1
    line_start = 0

    # Process text and yield tokens
    for mo in re.finditer(tok_regex, text):
        kind = mo.lastgroup
        value = mo.group()
        column = mo.start() - line_start
        if kind == "CHROMOSOME":
            value = value
        elif kind == "NUMBER":
            value = int(value)
        elif kind == "TAB":
            value = "\t"
        elif kind == "NEWLINE":
            line_start = mo.end()
            line_num += 1
        elif kind == "comment":
            line_start = mo.end()
            line_num += 1
        elif kind == "ANY":
            pass
        yield Token(kind, value, line_num, column)


class Record:
    """BED record"""

    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.region_id = f"{self.chrom}_{self.start}_{self.end}"

    def __str__(self):
        return f"{self.chrom}\t{str(self.start)}\t{str(self.end)}"


class ParserError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


def parse_BED(file_object):
    """This parser will take a file object as input, read the file
    and try to parse it as a bed file. If the content of the file
    doesn't adhere to the bed specification an Exception will be raised
    """
    record_values = cycle(["chrom", "start", "end"])
    sequence = cycle(["CHROMOSOME", "TAB", "NUMBER", "TAB", "NUMBER", "NEWLINE"])

    chrom = None
    start = None
    end = None
    next_sequence = next(sequence)
    next_record = next(record_values)

    for token in tokenize(file_object.read()):
        if next_sequence == "NEWLINE" and token.type != "NEWLINE":
            continue
        elif next_sequence == "CHROMOSOME" and token.type == "COMMENT":
            continue
        elif next_sequence == token.type:
            if next_record == "chrom" and token.type == "CHROMOSOME":
                chrom = token.value
                next_record = next(record_values)
            elif next_record == "start" and token.type == "NUMBER":
                start = token.value
                next_record = next(record_values)
            elif next_record == "end" and token.type == "NUMBER":
                end = token.value
                next_record = next(record_values)
                yield Record(chrom, start, end)
        else:
            logging.error("Error occured while parsing the bed file")
            raise ParserError(f"BED file format did not match for token {token}")
        next_sequence = next(sequence)


###########################################
# ---------| Output Writting |------------#
###########################################


def write_header(args, bam_report):
    header = "\n".join(
        [
            f"##fileDate={datetime.now().strftime('%Y%m%d')}",
            f"##source={__file__}",
            f"##bamFile={args.bam_file}",
            f"##bedFile={args.bed_file.name}",
            f"##minLength={args.min_length}",
            f"##maxLength={args.max_length}",
            f"##minMappingQuality={args.mapq}",
            f"##readsFetchedFromBam={bam_report.fetched_reads}",
            f"##readsPassedQualCheck={bam_report.reads_passed_qual_check}",
            f"##pairedReads={bam_report.paired_reads}",
            f"##pairedReadsPassedQualCheck={bam_report.paired_reads_passed_qual_check}",
            f"##pairedReadsEmitted={bam_report.paired_reads_yielded}",
            f"#{' '.join(['region'] + list(map(str, range(args.min_length, args.max_length+1))))}",
        ]
    )
    print(header)  # noqa: T001


def write_output(args, bam_report, length_dists, region_ids):
    write_header(args, bam_report)
    for region_id, length_dist in zip(region_ids, length_dists):
        print(" ".join([region_id] + list(map(str, length_dist))))  # noqa: T001


###########################################
# ------| Where The Magic Happens |-------#
###########################################


class ReadPair(NamedTuple):
    ref_name: str
    start: int
    end: int
    start_is_first: str
    length: int


class Report:
    def __init__(self, filename):
        self.filename = filename
        self.fetched_reads = 0
        self.reads_passed_qual_check = 0
        self.paired_reads = 0
        self.paired_reads_passed_qual_check = 0
        self.paired_reads_yielded = 0


class BAM:
    def __init__(self, filename):
        # Check if index exists, if not create an index file
        index_filename = f"{filename}.bai"
        if not os.path.exists(index_filename):
            logging.warning(f"No index file found ({index_filename}), generating...")
            pysam.index(filename)
        # Set class variables
        self.bam_file = pysam.AlignmentFile(filename, "rb")
        self.report = Report(filename)

    def _contains_clipping(self, read):
        return any(map(lambda c: c == "S" or c == "H", read.cigarstring))

    def fetch_readpairs(self, chrom, region_start, region_end, mapq):
        # Define variables
        unpaired_reads = {}

        # Fetch reads in region
        for read in self.bam_file.fetch(
            contig=chrom, start=region_start, stop=region_end
        ):
            self.report.fetched_reads += 1
            # Test read quality
            if (
                read.is_duplicate
                or read.is_secondary
                or read.is_supplementary
                or read.mapping_quality < mapq
                or self._contains_clipping(read)
            ):
                continue

            # Find query name
            query_name = read.query_name
            self.report.reads_passed_qual_check += 1

            if query_name not in unpaired_reads:
                # Store read until partner is found
                unpaired_reads[query_name] = (
                    read.reference_start,
                    read.reference_end,
                    read.is_reverse,
                    read.is_read1,
                )
            else:
                # fetch read partner
                mem_start, mem_end, mem_reverse, mem_is_read1 = unpaired_reads[
                    query_name
                ]
                del unpaired_reads[query_name]
                self.report.paired_reads += 2
                # Test readpair quality
                if (
                    mem_start is None
                    or mem_end is None
                    or read.reference_start is None
                    or read.reference_end is None
                    or read.is_reverse == mem_reverse
                ):
                    continue
                self.report.paired_reads_passed_qual_check += 2
                # find fragment ends and first mate
                if read.is_reverse:
                    start = mem_start
                    end = read.reference_end
                    start_is_first = mem_is_read1
                else:
                    start = read.reference_start
                    end = mem_end
                    start_is_first = not mem_is_read1
                # Extract read length
                length = abs(read.template_length)

                # Last quality check
                if start >= end or length == 0:
                    continue
                if length != end - start:
                    logging.warning(
                        f"Aligner specified length and length computed from fragmen endpoints does not match. Template length: {length}, Computed length:{end-start}."
                    )

                self.report.paired_reads_yielded += 2

                yield ReadPair(
                    read.reference_name, int(start), int(end), start_is_first, length
                )


def compute_length_distributions(args):
    # Parse bed file
    logging.debug("Parsing bed file...")
    regions = list(parse_BED(args.bed_file))
    n_regions = len(regions)
    logging.debug(f"Bed file parsed found {n_regions} regions")

    # Prepare bam file
    bam = BAM(args.bam_file)

    # Start fragment extraction
    logging.info("Start extracting fragments from bam file...")
    length_dists = np.zeros(
        (n_regions, args.max_length - args.min_length + 1), dtype=np.uint
    )
    region_ids = list()
    for i, region in enumerate(regions):
        log_progress(i, n_regions)
        region_ids.append(region.region_id)

        for readpair in bam.fetch_readpairs(
            region.chrom, region.start, region.end, args.mapq
        ):
            length = readpair.length
            if args.min_length <= length <= args.max_length:
                length_dists[i, length - args.min_length] += 1

    logging.info("Finished extration, start writting output")
    write_output(args, bam.report, length_dists, region_ids)
    logging.info("Finished writting output")


###########################################
# ----| Argument Parsing And Logging|-----#
###########################################


def argument_parsing():
    parser = argparse.ArgumentParser(
        description="A script for extracting the length distribution of ctDNA WGS data"
    )
    parser.add_argument(
        "bam_file", help="File path to the bam file that should be processed"
    )
    parser.add_argument(
        "bed_file",
        type=argparse.FileType("r"),
        help="File path to the BED file defining the genomic regions that should be processed",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=500,
        help="Maximum fragment length to include",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=1,
        help="Minimum fragment length to include",
    )
    parser.add_argument(
        "--mapq",
        type=int,
        default=20,
        help="Minimum mapping quality",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="File path to where the result should be outputtet. Otherwise stdout is used",
    )
    parser.add_argument(
        "-v",
        action="store_true",
        help="Set logging level to debug",
    )
    parser.add_argument(
        "--no-cols",
        action="store_true",
        help="Remove color formatting",
    )
    return parser.parse_args()


def setup_logging(args):
    level = logging.DEBUG if args.v else logging.INFO
    GREEN = "" if args.no_cols else "\033[92m"
    WHITE = "" if args.no_cols else "\033[00m"
    logging.basicConfig(
        format=f"{GREEN}[%(levelname)s: %(asctime)s]:{WHITE} %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=level,
    )


def log_progress(i, N):
    percentage_steps = 10
    if N < percentage_steps:
        return
    if i % (N // percentage_steps) == 0:
        percentage = (i / N) * 100
        percentage_rounded = int(
            math.ceil(percentage / percentage_steps) * percentage_steps
        )
        logging.info("Extracting fragments... {}%".format(percentage_rounded))


###########################################
# ---------| Script Entrypoint |----------#
###########################################

BANNER = b"CiAgICAgICAgICAgICBfXyAgICAgIF9fX19fX18gICBfXyAgICBfXyAgIF9fX19fXyAgICAgICAgICAgX18gICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBfXyAgICAgIF9fICAgICAgIAogICAgICAgICAgICB8ICBcICAgIHwgICAgICAgXCB8ICBcICB8ICBcIC8gICAgICBcICAgICAgICAgfCAgXCAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB8ICBcICAgIHwgIFwgICAgICAKICBfX19fX19fICBffCAkJF8gICB8ICQkJCQkJCRcfCAkJFwgfCAkJHwgICQkJCQkJFwgICAgICAgIHwgJCQgIF9fX19fXyAgIF9fX19fX18gICAgX19fX19fICBffCAkJF8gICB8ICQkX19fXyAgCiAvICAgICAgIFx8ICAgJCQgXCAgfCAkJCAgfCAkJHwgJCQkXHwgJCR8ICQkX198ICQkIF9fX19fXyB8ICQkIC8gICAgICBcIHwgICAgICAgXCAgLyAgICAgIFx8ICAgJCQgXCAgfCAkJCAgICBcIAp8ICAkJCQkJCQkIFwkJCQkJCQgIHwgJCQgIHwgJCR8ICQkJCRcICQkfCAkJCAgICAkJHwgICAgICBcfCAkJHwgICQkJCQkJFx8ICQkJCQkJCRcfCAgJCQkJCQkXCAkJCQkJCQgIHwgJCQkJCQkJCQKfCAkJCAgICAgICAgfCAkJCBfXyB8ICQkICB8ICQkfCAkJFwkJCAkJHwgJCQkJCQkJCQgXCQkJCQkJHwgJCR8ICQkICAgICQkfCAkJCAgfCAkJHwgJCQgIHwgJCQgfCAkJCBfXyB8ICQkICB8ICQkCnwgJCRfX19fXyAgIHwgJCR8ICBcfCAkJF9fLyAkJHwgJCQgXCQkJCR8ICQkICB8ICQkICAgICAgICB8ICQkfCAkJCQkJCQkJHwgJCQgIHwgJCR8ICQkX198ICQkIHwgJCR8ICBcfCAkJCAgfCAkJAogXCQkICAgICBcICAgXCQkICAkJHwgJCQgICAgJCR8ICQkICBcJCQkfCAkJCAgfCAkJCAgICAgICAgfCAkJCBcJCQgICAgIFx8ICQkICB8ICQkIFwkJCAgICAkJCAgXCQkICAkJHwgJCQgIHwgJCQKICBcJCQkJCQkJCAgICBcJCQkJCAgXCQkJCQkJCQgIFwkJCAgIFwkJCBcJCQgICBcJCQgICAgICAgICBcJCQgIFwkJCQkJCQkIFwkJCAgIFwkJCBfXCQkJCQkJCQgICBcJCQkJCAgXCQkICAgXCQkCiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICB8ICBcX198ICQkICAgICAgICAgICAgICAgICAgIAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIFwkJCAgICAkJCAgICAgICAgICAgICAgICAgICAKICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgXCQkJCQkJCAgICAgICAgICAgICAgICAgICAgCg=="
BANNER = base64.b64decode(BANNER).decode()


def main():
    args = argument_parsing()
    setup_logging(args)
    print(BANNER, file=sys.stderr)  # noqa: T001
    compute_length_distributions(args)


if __name__ == "__main__":
    main()
