#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
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
        "--outdir", dest="output_dir", type=Path, required=True, help="Output directory"
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


def get_coverage(hit_lf: pl.LazyFrame, target: str):
    # takes the sub dataframe corresponding to this target
    df = hit_lf.filter(pl.col("sseqid") == target).collect()

    if df.is_empty():
        return pd.DataFrame()

    if len(df.select("slen").unique()) > 1:
        raise ValueError("Multiple lengths")

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
    output_dir: Path,
    log: bool = False,
):
    fig, ax = plt.subplots(figsize=FIGZISE)

    common_plot_params = dict(
        kind="area", ax=ax, alpha=0.6, linewidth=0, stacked=False, xlabel=target
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

    # cleaning target name
    dirname = output_dir / target.replace("/", "_")
    Path(dirname).mkdir(parents=True)
    outfile = f"{dirname}/log.png" if log else f"{dirname}/raw.png"

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

    hit_lf = pl.scan_csv(
        args.blast_hit_file, separator="\t", schema=BLAST_OUTPUT_COL_SCHEMA
    )
    chimera_lf = pl.scan_csv(args.chimera_file)

    # dividing hit table into chimera and non-chimera hits
    chimera_qseqids = chimera_lf.select("qseqid").collect().to_series().to_list()
    chimera_hit_lf = hit_lf.filter(pl.col("qseqid").is_in(chimera_qseqids))
    non_chimera_hit_lf = hit_lf.filter(~pl.col("qseqid").is_in(chimera_qseqids))

    # looping through target sequences and making coverage plot
    target_sequences = hit_lf.select("sseqid").collect().unique().to_series().to_list()

    non_chimera_coverage = {}
    chimera_coverage = {}
    logger.info("Computing coverages for each target sequence")
    for target in tqdm(target_sequences):
        coverage_df = get_coverage(chimera_hit_lf, target)
        # we don't want to keep targets for which there is not chimera
        if coverage_df.empty:
            continue
        chimera_coverage[target] = coverage_df
        non_chimera_coverage[target] = get_coverage(non_chimera_hit_lf, target)

    # plotting images
    logger.info("Making plots")
    for target in tqdm(chimera_coverage):
        plot_coverage(
            non_chimera_coverage[target],
            chimera_coverage[target],
            target,
            args.output_dir,
        )
        plot_coverage(
            non_chimera_coverage[target],
            chimera_coverage[target],
            target,
            args.output_dir,
            log=True,
        )

    logger.info("Done")
