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
    parser = argparse.ArgumentParser(
        description="Compute general statistics from count data for each sample"
    )
    parser.add_argument(
        "--hits", type=Path, dest="blast_hits_file", required=True
    )
    parser.add_argument(
        "--out", type=Path, dest="outfile", required=True
    )
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()
    logger.info(f"Parsing {args.blast_hits_file}")

    df = pd.read_csv(args.blast_hits_file)

    family_to_nb_hits = df.groupby("family")["nb_hits"].agg(list).to_dict()

    # add None values when necessary
    max_len = max(len(nb_hits_list) for nb_hits_list in family_to_nb_hits.values())

    family_to_nb_hits_aligned = {
        family: nb_hits_list + [None] * (max_len - len(nb_hits_list))
        for family, nb_hits_list in family_to_nb_hits.items()
    }

    formatted_df = pd.DataFrame(family_to_nb_hits_aligned)

    # in case families are all taxids, cas them to string
    with open(args.outfile, "w", newline="") as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        writer.writerow(formatted_df.columns)

        formatted_df.to_csv(f, index=False)
        logger.info("Done")
