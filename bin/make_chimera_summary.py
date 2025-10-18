#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import pandas as pd
import argparse
from pathlib import Path
import logging

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

GROUPBY_COLS = {"species": "species.taxid", "family": "family"}

COLS_TO_REMOVE = [
    "species.taxid",
    "family",
    "qstart_1",
    "qend_1",
    "sstart_1",
    "send_1",
    "qstart_2",
    "qend_2",
    "sstart_2",
    "send_2",
    "coordinate_in_1",
    "coordinate_in_2",
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=Path, dest="chimeras_file", required=True)
    parser.add_argument("--groupby", type=str, dest="col_for_grouping", required=True)
    parser.add_argument("--out", type=Path, dest="outfile", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    chimeras_df = pd.read_csv(args.chimeras_file)
    col_for_grouping = GROUPBY_COLS[args.col_for_grouping]

    logger.info("Converting all values to numeric when possible")
    # force convert to numeric whenever possible
    chimeras_df = chimeras_df.apply(pd.to_numeric, errors="ignore")

    logger.info("CSelecting all numerical columns")
    # selecting all float, int and bool columns
    cols_to_keep = [col_for_grouping] + [
        col
        for col in chimeras_df.columns
        if chimeras_df[col].dtype != "object" and col not in COLS_TO_REMOVE
    ]
    numeric_df = chimeras_df[cols_to_keep]

    logger.info("Applying mean")
    # grouping by family or taxid & computing mean
    means_df = numeric_df.groupby(col_for_grouping).mean().reset_index()
    print(means_df["slen_2"])
    logger.info("Exporting")
    means_df.to_csv(args.outfile, index=False, header=True)

    logger.info("Done")
