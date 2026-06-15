#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

import polars as pl

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)


CHIMERA_COL_SCHEMA = pl.Schema(
    {
        "qseqid": pl.String(),
        "sseqid_1": pl.String(),
        "pident_1": pl.Float64(),
        "length_1": pl.UInt32(),
        "mismatch_1": pl.UInt32(),
        "gapopen_1": pl.UInt32(),
        "qstart_1": pl.UInt32(),
        "qend_1": pl.UInt32(),
        "sstart_1": pl.UInt64(),
        "send_1": pl.UInt64(),
        "qlen_1": pl.UInt32(),
        "slen_1": pl.UInt64(),
        "evalue_1": pl.Float64(),
        "bitscore_1": pl.Float64(),
        "sseqid_2": pl.String(),
        "pident_2": pl.Float64(),
        "length_2": pl.UInt32(),
        "mismatch_2": pl.UInt32(),
        "gapopen_2": pl.UInt32(),
        "qstart_2": pl.UInt32(),
        "qend_2": pl.UInt32(),
        "sstart_2": pl.UInt64(),
        "send_2": pl.UInt64(),
        "qlen_2": pl.UInt32(),
        "slen_2": pl.UInt64(),
        "evalue_2": pl.Float64(),
        "bitscore_2": pl.Float64(),
        "total_coverage": pl.Int32(),
        "overlap_length": pl.Int32(),
        "coverage_1_only": pl.Int32(),
        "coverage_2_only": pl.Int32(),
        "chimeric": pl.Boolean(),
        "coordinate_in_1": pl.Float64(),
        "coordinate_in_2": pl.Float64(),
        "family": pl.String(),
        "species.taxid": pl.UInt64(),
        "srr": pl.String(),
    }
)

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get chimera and non-chimera read coverage for each target"
    )
    parser.add_argument(
        "--chimeras",
        type=Path,
        dest="chimeras_files",
        nargs="+",
        required=True,
        help="Chimera table file",
    )
    parser.add_argument(
        "--out",
        type=Path,
        dest="out_file",
        required=True,
        help="Output file path",
    )
    return parser.parse_args()


def parse_chimeras_file(file: Path) -> pl.LazyFrame:
    return pl.scan_csv(file, has_header=True, schema=CHIMERA_COL_SCHEMA)


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    # creating output file anyway
    Path(args.out_file).touch()

    with open(args.out_file, "a") as f:
        for i, file in enumerate(args.chimeras_files):
            chimeras_lf = parse_chimeras_file(file)
            try:
                if i == 0:
                    chimeras_lf.sink_csv(f, include_header=True)
                else:
                    chimeras_lf.sink_csv(f, include_header=False)
            except pl.exceptions.NoDataError:
                logger.info(f"No data in file {file}")


    logger.info("Done")
