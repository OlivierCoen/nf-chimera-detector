#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import pandas as pd
import argparse
from pathlib import Path
import logging

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--files", type=str, dest="chimeras_files", required=True
    )
    parser.add_argument(
        "--out", type=Path, dest="outfile", required=True
    )
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()
    chimeras_dfs = [pd.read_csv(file) for file in args.chimeras_files.split(" ")]
    chimeras_df = pd.concat(chimeras_dfs, ignore_index=True)

    chimeras_df.to_csv(args.outfile, index=False, header=True)
    print(chimeras_df)
    logger.info("Done")

