#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path
from tqdm import tqdm
import polars as pl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

BLAST_OUTPUT_COL_SCHEMA = pl.Schema(
    {
        "qseqid": pl.String(),
        "sseqid": pl.String(),
        "pident": pl.Float64(),
        "length": pl.UInt32(),
        "mismatch": pl.UInt32(),
        "gapopen": pl.UInt32(),
        "qstart": pl.UInt32(),
        "qend": pl.UInt32(),
        "sstart": pl.UInt32(),
        "send": pl.UInt32(),
        "qlen": pl.UInt64(),
        "slen": pl.UInt64(),
        "evalue": pl.Float64(),
        "bitscore": pl.Float64(),
    }
)

CHIMERA_COL_SCHEMA = {
    "qseqid": "string",
    "sseqid_1": "string",
    "pident_1": "Float64",
    "length_1": "Int64",
    "mismatch_1": "Int64",
    "gapopen_1": "Int64",
    "qstart_1": "Int64",
    "qend_1": "Int64",
    "sstart_1": "Int64",
    "send_1": "Int64",
    "qlen_1": "Int64",
    "slen_1": "Int64",
    "evalue_1": "Float64",
    "bitscore_1": "Float64",
    "sseqid_2": "string",
    "pident_2": "Float64",
    "length_2": "Int64",
    "mismatch_2": "Int64",
    "gapopen_2": "Int64",
    "qstart_2": "Int64",
    "qend_2": "Int64",
    "sstart_2": "Int64",
    "send_2": "Int64",
    "qlen_2": "Int64",
    "slen_2": "Int64",
    "evalue_2": "Float64",
    "bitscore_2": "Float64",
    "total_coverage": "Int64",
    "overlap_length": "Int64",
    "coverage_1_only": "Int64",
    "coverage_2_only": "Int64",
    "chimeric": "boolean",
    "coordinate_in_1": "Int64",
    "coordinate_in_2": "Int64",
    "family": "string",
    "species.taxid": "Int64",
    "srr": "string",
}

FIGZISE = (12, 6)

CHIMERA_COLOR = "green"
NON_CHIMERA_COLOR = "grey"

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
        "--chimeras",
        dest="chimera_file",
        type=Path,
        required=True,
        help="File containing chimera table",
    )
    parser.add_argument(
        "--srr",
        dest="srr_id",
        type=str,
        required=True,
        help="SRR ID",
    )
    return parser.parse_args()


def get_occurences(df: pl.DataFrame):
    target_length = df.select("slen").unique().to_series().item(0)
    fields = [str(i) for i in range(target_length)]

    return (
        df.with_columns(
            nb_left_uncovered=(
                pl.col("sstart") - 1
            ),  # nb of bases strictly before sstart
            nb_covered=(
                pl.col("send") - pl.col("sstart") + 1
            ),  # nb of bases between ssart and send
            nb_right_uncovered=(
                pl.lit(target_length) - pl.col("send")
            ),  # nb of bases stricly after send
        )
        .with_columns(
            left_uncovered=pl.int_ranges("nb_left_uncovered").list.eval(
                pl.element() * 0
            ),
            # list of zeros of length nb_left_uncovered
            covered=pl.int_ranges("nb_covered").list.eval(
                (pl.element() >= 0).cast(pl.Int64)
            ),
            # trick to have list of ones of length nb_covered
            right_uncovered=pl.int_ranges("nb_right_uncovered").list.eval(
                pl.element() * 0
            ),
            # list of zeros of length nb_right_uncovered
        )
        .with_columns(
            occurences=pl.concat_list(["left_uncovered", "covered", "right_uncovered"])
            # concatenation of the threee lists: list of zeros / ones of length target_length
        )
        .select(
            pl.col("occurences").list.to_struct(
                fields=fields
            )  # selecting only this column and preparing for unnesting
        )
        .unnest("occurences")
        # turning "occurences" column (column of lists) into a dataframe: each column correspond to a position on the subject sequence
    )


def swap_start_end_if_necessary(df: pl.DataFrame):
    cond = pl.col("sstart") > pl.col("send")
    return df.with_columns(
        [
            pl.when(cond)
            .then(pl.col("send"))
            .otherwise(pl.col("sstart"))
            .alias("sstart"),
            pl.when(cond)
            .then(pl.col("sstart"))
            .otherwise(pl.col("send"))
            .alias("send"),
        ]
    )


def get_coverage(hit_df: pl.DataFrame, target: str):
    # takes the sub dataframe corresponding to this target
    df = hit_df.filter(pl.col("sseqid") == target)

    if df.is_empty():
        return pd.DataFrame()

    # if multiple subject lengths are found, log a warning
    if len(df.select("slen").unique()) > 1:
        logger.warning(f"Multiple subject lengths for target {target}:\n {df}")

    # swapping start and end if on the reverse strand
    df = swap_start_end_if_necessary(df)

    # compute a dataframe of shape (nb of hits * subject length) with 0 and 1
    # 1 means that a hit was included at that position, and 0 that it was not at this position
    occurence_df = get_occurences(df)

    # returns the sum of all occurences for each position
    return (
        occurence_df.select(
            pl.sum("*")
        )  # for each column (each position), suming all occurences (1 or 0)
        .transpose()
        .to_series()
        .to_pandas()
    )


def plot_coverage(
    non_chimera_coverage: pd.Series,
    chimera_coverage: pd.Series,
    target: str,
    srr_id: str,
    log: bool = False,
):
    fig, ax = plt.subplots(figsize=FIGZISE)

    xlabel = f"SRR: {srr_id} / TE: {target} "
    common_plot_params = dict(
        kind="area", ax=ax, alpha=0.6, linewidth=0, stacked=False, xlabel=xlabel
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PLOTTING NON CHIMERAS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if not non_chimera_coverage.empty:
        s1 = (
            non_chimera_coverage.apply(lambda x: 10 * np.log10(x + 1))
            if log
            else non_chimera_coverage
        )
        s1.plot(label="Non chimeras", color=NON_CHIMERA_COLOR, **common_plot_params)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PLOTTING CHIMERAS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # chimera_coverage is necessarily not empty because we filtered it beforehand
    s2 = (
        chimera_coverage.apply(lambda x: 10 * np.log10(x + 1))
        if log
        else chimera_coverage
    )
    s2.plot(label="Chimeras", color=CHIMERA_COLOR, **common_plot_params)

    cleaned_target = target.replace("/", "_")
    outfile = (
        f"{srr_id}.{cleaned_target}.log.png"
        if log
        else f"{srr_id}.{cleaned_target}.raw.png"
    )

    plt.savefig(outfile, bbox_inches="tight")
    # closes current figure window to save memory
    plt.close()


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    hit_lf = pl.read_csv(
        args.blast_hit_file,
        separator="\t",
        schema=BLAST_OUTPUT_COL_SCHEMA,
        has_header=False,
    )

    try:
        chimera_df = pd.read_csv(
            args.chimera_file,
            dtype=CHIMERA_COL_SCHEMA,
        )
    except pd.errors.EmptyDataError as e:
        logger.warning(f"No data in chimera file: {e}")
        sys.exit(0)

    # dividing hit table into chimera and non-chimera hits
    chimera_qseqids = chimera_df["qseqid"].unique().tolist()
    chimera_hit_df = hit_lf.filter(pl.col("qseqid").is_in(chimera_qseqids))
    non_chimera_hit_df = hit_lf.filter(~pl.col("qseqid").is_in(chimera_qseqids))

    # looping through target sequences and making coverage plot
    target_sequences = hit_lf.select("sseqid").unique().to_series().to_list()

    non_chimera_coverage = {}
    chimera_coverage = {}
    logger.info("Computing coverages for each target sequence")
    for target in tqdm(target_sequences):
        coverage_df = get_coverage(chimera_hit_df, target)
        # we don't want to keep targets for which there is not chimera
        if coverage_df.empty:
            continue
        chimera_coverage[target] = coverage_df
        non_chimera_coverage[target] = get_coverage(non_chimera_hit_df, target)

    # plotting images
    logger.info("Making plots")
    for target in tqdm(chimera_coverage):
        plot_coverage(
            non_chimera_coverage[target],
            chimera_coverage[target],
            target,
            args.srr_id,
        )
        plot_coverage(
            non_chimera_coverage[target],
            chimera_coverage[target],
            target,
            args.srr_id,
            log=True,
        )

    logger.info("Done")
