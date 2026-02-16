#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
from tqdm import tqdm

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

CHUNKSIZE = 100000
MAX_NUMBER_HITS_CONSIDERED = 100000

CHIMERA_COLOR = "green"
NON_CHIMERA_COLOR = "grey"

COMMON_PLOT_PARAMS = dict(kind="area", alpha=0.6, linewidth=0, stacked=False)

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

    """
    with open("log.txt", "a") as f:
        msg = f"target: {df['sseqid'].item(0)} || nb lines: {len(df)}"
        f.write(msg)
    """

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


def get_chunk_coverage(df: pl.DataFrame) -> pl.DataFrame:
    # swapping start and end if on the reverse strand
    df = swap_start_end_if_necessary(df)

    # compute a dataframe of shape (nb of hits * subject length) with 0 and 1
    # 1 means that a hit was included at that position, and 0 that it was not at this position
    occurence_df = get_occurences(df)

    # returns the sum of all occurences for each position
    return occurence_df.select(
        pl.sum("*")
    ).transpose()  # for each column (each position), suming all occurences (1 or 0)


def get_coverage(
    hit_lf: pl.LazyFrame, target: str, is_chimera: bool
) -> tuple[pl.DataFrame, bool]:
    # takes the sub dataframe corresponding to this target
    lf = hit_lf.filter(pl.col("sseqid") == target)

    if lf.limit(1).collect().is_empty():
        return pl.DataFrame(), False

    # if multiple subject lengths are found, log a warning
    if len(lf.select("slen").collect().unique()) > 1:
        logger.warning(f"Multiple subject lengths for target {target}")

    # get the total number of rows
    n_rows = lf.select(pl.len()).collect().item()

    is_subsampled = False
    if not is_chimera and n_rows > MAX_NUMBER_HITS_CONSIDERED:
        logger.warning(
            f"Too many hits for target {target}: {n_rows}. Subsampling to {MAX_NUMBER_HITS_CONSIDERED}"
        )
        indices = sorted(
            np.random.choice(n_rows, size=MAX_NUMBER_HITS_CONSIDERED, replace=False)
        )
        lf = (
            lf.with_row_index()
            .filter(pl.col("index").is_in(indices))
            .select(pl.exclude("index"))
        )
        is_subsampled = True

    # perform coverage computation on chunks to limit memory usage
    coverage_dfs = []
    for i in range(0, n_rows, CHUNKSIZE):
        chunk = lf.slice(i, CHUNKSIZE).collect()
        if chunk.is_empty():
            break
        coverage_dfs.append(get_chunk_coverage(chunk))

    df = pl.concat(coverage_dfs, how="vertical")

    if is_chimera:
        df = df.rename({"column_0": "chimera"})
    else:
        df = df.rename({"column_0": "non chimera"})

    return df, is_subsampled


def get_series_to_plot(df: pl.DataFrame, log: bool) -> pd.Series:
    s = df[df.columns[0]].to_pandas()
    return s.apply(lambda x: 10 * np.log10(x + 1)) if log else s


def get_outfile(target: str, srr_id: str, log: bool, is_subsampled: bool) -> str:
    cleaned_target = target.replace("/", "_")
    subsampling_status = ".subsampled" if is_subsampled else ""
    return (
        f"{srr_id}.{cleaned_target}{subsampling_status}.log.png"
        if log
        else f"{srr_id}.{cleaned_target}{subsampling_status}.raw.png"
    )


def plot_coverage(
    non_chimera_coverage_df: pl.DataFrame,
    chimera_coverage_df: pl.DataFrame,
    target: str,
    srr_id: str,
    is_subsampled: bool,
    log: bool = False,
):
    fig, ax = plt.subplots(figsize=FIGZISE)

    xlabel = f"SRR: {srr_id} / TE: {target} "

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PLOTTING NON CHIMERAS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if not non_chimera_coverage_df.is_empty():
        non_chimera_coverage = get_series_to_plot(non_chimera_coverage_df, log)
        non_chimera_coverage.plot(
            y="non chimera",
            label="Non chimeras",
            color=NON_CHIMERA_COLOR,
            ax=ax,
            xlabel=xlabel,
            **COMMON_PLOT_PARAMS,
        )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PLOTTING CHIMERAS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # chimera_coverage is necessarily not empty because we filtered it beforehand
    chimera_coverage = get_series_to_plot(chimera_coverage_df, log)
    chimera_coverage.plot(
        label="Chimeras",
        color=CHIMERA_COLOR,
        ax=ax,
        xlabel=xlabel,
        **COMMON_PLOT_PARAMS,
    )

    outfile = get_outfile(target, srr_id, log, is_subsampled)

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

    logger.info("Dividing hit table into chimera and non-chimera hits")
    chimera_qseqids = chimera_df["qseqid"].unique().tolist()
    chimera_hit_lf = hit_lf.filter(pl.col("qseqid").is_in(chimera_qseqids))
    non_chimera_hit_lf = hit_lf.filter(~pl.col("qseqid").is_in(chimera_qseqids))

    # looping through target sequences and making coverage plot
    logger.info("Extracting target sequences")
    target_sequences = hit_lf.select("sseqid").unique().collect().to_series().to_list()

    logger.info("Computing coverages for each target sequence and making plots")
    for target in tqdm(target_sequences):
        # getting coverage for chimera reads
        chimera_coverage_df, is_subsampled = get_coverage(
            chimera_hit_lf, target, is_chimera=True
        )

        # we don't want to keep targets for which there is not chimera
        if chimera_coverage_df.is_empty():
            continue

        non_chimera_coverage_df, is_subsampled = get_coverage(
            non_chimera_hit_lf, target, is_chimera=False
        )

        plot_coverage(
            non_chimera_coverage_df,
            chimera_coverage_df,
            target,
            args.srr_id,
            is_subsampled,
        )

        plot_coverage(
            non_chimera_coverage_df,
            chimera_coverage_df,
            target,
            args.srr_id,
            is_subsampled,
            log=True,
        )

    logger.info("Done")
