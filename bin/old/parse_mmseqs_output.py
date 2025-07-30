#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import pandas as pd
from pathlib import Path

import logging

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

HEADER = ["query", "target", "pident", "alilen", "mismatches", "gapopens", "qstart", "qend", "tstart", "tend", "evalue", "bits"]

OUTFILE = "status.txt"

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse and filter MMSEQS output file"
    )
    parser.add_argument(
        "--file", type=Path, dest="mmseqs_output_file", required=True
    )
    parser.add_argument(
        "--max-evalue", type=float, dest="max_evalue", required=True
    )
    return parser.parse_args()

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()
    mmseqs_output_file = args.mmseqs_output_file

    logger.info(f"Parsing {mmseqs_output_file}")
    df = pd.read_csv(mmseqs_output_file, names=HEADER, sep='\t')

    with open(OUTFILE, 'w') as fout:
        if (df["evalue"] <= args.max_evalue).any():
            logger.info("Found results with evalue lower or equal than max threshold")
            fout.write("PASS")
        else:
            logger.info("All parsed evalues are higher than max threshold")
            fout.write("FAIL")

    logger.info("Done")
