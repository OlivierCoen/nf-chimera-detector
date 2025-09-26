#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import pandas as pd
import argparse
from pathlib import Path
import csv
import logging

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data", type=Path, dest="data_file", required=True
    )
    parser.add_argument(
        "--out", type=Path, dest="outfile", required=True
    )
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()
    logger.info(f"Parsing {args.data_file}")

    df = pd.read_csv(args.data_file, sep='\t')

    family_to_data = df.groupby("family")["data"].agg(list).to_dict()

    # add None values when necessary
    max_len = max(len(nb_hits_list) for nb_hits_list in family_to_data.values())

    family_to_data_aligned = {
        family: nb_hits_list + [None] * (max_len - len(nb_hits_list))
        for family, nb_hits_list in family_to_data.items()
    }

    formatted_df = pd.DataFrame(family_to_data_aligned)

    with open(args.outfile, "w", newline="") as f:
        # print header separately with quotes
        writer = csv.writer(f, quoting=csv.QUOTE_ALL, delimiter ='\t')
        writer.writerow(formatted_df.columns)

        formatted_df.to_csv(f, index=False, header=False, sep='\t')
        logger.info("Done")
