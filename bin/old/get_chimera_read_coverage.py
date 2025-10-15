#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path
from tqdm import tqdm
import pandas as pd
from functools import partial
import matplotlib.pyplot as plt


logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

BLAST_OUTPUT_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "qlen", "slen", "evalue", "bitscore"]
FIGZISE = (12, 6)

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
        "--hits", type=Path, dest="blast_hit_file", required=True, help="Blast hit file"
    )
    parser.add_argument(
        "--chimeras", dest="chimera_file", type=Path, required=True, help="File containing chimera table"
    )
    return parser.parse_args()


def fill_positions(start: int, end: int, te_length: int) -> list[int]:
    if start > end:
        start, end = end, start
    left_uncovered_positions = [0] * (start - 1)
    covered_positions = [1] * (end - start + 1)
    right_uncovered_positions = [0] * (te_length - end)
    return left_uncovered_positions + covered_positions + right_uncovered_positions


def get_coverage(hit_df: pd.DataFrame, target: str) -> pd.Series:
    te_hit_df = hit_df[hit_df.sseqid == target]
    if te_hit_df.empty:
        return None
    if len(te_hit_df.slen.unique()) > 1:
        raise ValueError("Multiple lengths")
    te_length = te_hit_df.slen.unique()[0]
    func = partial(fill_positions, te_length=te_length)
    occurence_df =  te_hit_df.apply(lambda row: func(row.sstart, row.send), axis=1, result_type="expand")
    return occurence_df.sum()


def plot_coverage(
        non_chimera_coverage: pd.Series,
        chimera_coverage: pd.Series,
        outfile: str
    ):
    fig, ax = plt.subplots(figsize=FIGZISE)
    non_chimera_coverage.plot(kind='area', ax=ax, alpha=0.4, label='Non chimeras', linewidth=0, color="grey")
    chimera_coverage.plot(kind='area', ax=ax, alpha=0.4, label='Chimera', linewidth=0, color="red")
    plt.savefig(outfile)
    # closes current figure window to save memory
    plt.close()


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    hit_df = pd.read_csv(args.blast_hit_file, sep='\t', names=BLAST_OUTPUT_COLS)
    chimera_df = pd.read_csv(args.chimera_file)

    # dividing hit table into chimera and non-chimera hits
    chimera_qseqids = chimera_df.qseqid.tolist()
    chimera_hit_df = hit_df[hit_df.qseqid.isin(chimera_qseqids)]
    non_chimera_hit_df = hit_df[~hit_df.qseqid.isin(chimera_qseqids)]
    te_sequences = ['Wasmannia_auropunctata:rnd-5_family-1144#LTR/Viper']
    # looping through target sequences and making coverage plot
    for target in tqdm(te_sequences):
        non_chimera_coverage = get_coverage(non_chimera_hit_df, target)
        chimera_coverage = get_coverage(chimera_hit_df, target)
        # plotting image
        outfile = f"{target.replace('/', '_')}.png"
        plot_coverage(non_chimera_coverage, chimera_coverage, outfile)

    logger.info("Done")
