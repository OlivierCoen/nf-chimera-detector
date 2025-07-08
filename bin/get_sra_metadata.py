#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import pandas as pd
import subprocess
import argparse
from tqdm import tqdm
import logging

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute general statistics from count data for each sample"
    )
    parser.add_argument(
        "--taxon-id", type=int, dest="taxid", required=True, help="Taxon ID on NCBI Taxonomy"
    )
    return parser.parse_args()


def fetch_sra_experiment_uids_from_taxid(taxid: int) -> list[str]:
    command = f'esearch -db sra -query "txid{taxid}" | efetch -format uid'
    result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
    return result.stdout.strip().split("\n")

def fetch_sra_experiment_details(sra_uid: str) -> dict[str, str]:
    if sra_uid == "":
        logger.warning("Empty SRA ID. Skipping")
        return {}

    command = f'esearch -db sra -query "{sra_uid}[uid]" | efetch'
    result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
    lines = result.stdout.split("\n")
    keys = lines[0].split(",")
    exp_data = lines[1].split(",")
    return {key: value for key, value in zip(keys, exp_data)}


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    logger.info(f"Fetching sra experiment UIDs from taxon ID {args.taxid}")
    sra_uids = fetch_sra_experiment_uids_from_taxid(args.taxid)

    logger.info("Fetching sra experiment metadata for each SRA experiment ID")
    exp_dicts = []
    for sra_uid in tqdm(sra_uids, total=len(sra_uids)):
        exp_dict = fetch_sra_experiment_details(sra_uid)
        exp_dicts.append(exp_dict)

    logger.info("Formatting data and exporting...")
    df = pd.DataFrame.from_dict(exp_dicts)
    srrs = df["Run"].tolist()

    srr_outfile = f"{args.taxid}_srrs.txt"
    with open(srr_outfile, "w") as fout:
        fout.write("\n".join(srrs))

    sra_metadata_outfile = f"{args.taxid}_sra_metadata.csv"
    df.to_csv(sra_metadata_outfile, index=False, header=True)

    logger.info("Done")
