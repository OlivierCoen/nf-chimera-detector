#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import pandas as pd
import argparse
from pathlib import Path
import logging

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

GROUPBY_COLS = {
    "species": "species.taxid",
    "family": "family"
}

COLS_TO_REMOVE = [
    'species.taxid',
    'family',
    'chimeric',
    'qStart',
    'qEnd',
    'sStart',
    'sEnd',
    'qStart.s',
    'qEnd.s',
    'sStart.s',
    'sEnd.s',
    'insertionCoord',
    'chimericPoint'
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--files", type=str, dest="chimeras_files", required=True
    )
    parser.add_argument(
        "--groupby", type=str, dest="col_for_grouping", required=True
    )
    parser.add_argument(
        "--out", type=Path, dest="outfile", required=True
    )
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()
    chimeras_dfs = [pd.read_csv(file) for file in args.chimeras_files.split(" ")]
    col_for_grouping = GROUPBY_COLS[args.col_for_grouping]

    chimeras_df = pd.concat(chimeras_dfs, ignore_index=True)

    # force convert to numeric whenever possible
    chimeras_df = chimeras_df.apply(pd.to_numeric, errors='ignore')

    # selecting all float, int and bool columns
    cols_to_keep = [col_for_grouping] + [
        col for col in chimeras_df.columns
        if chimeras_df[col].dtype != 'object'
        and col not in COLS_TO_REMOVE
    ]
    numeric_df = chimeras_df[cols_to_keep]

    # grouping by family or taxid & computing mean
    means_df = numeric_df.groupby(col_for_grouping).mean().reset_index()

    means_df.to_csv(args.outfile, index=False, header=True)

