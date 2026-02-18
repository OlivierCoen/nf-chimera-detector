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
        "sstart": pl.UInt64(),
        "send": pl.UInt64(),
        "qlen": pl.UInt32(),
        "slen": pl.UInt64(),
        "evalue": pl.Float64(),
        "bitscore": pl.Float64(),
    }
)

CHIMERA_COL_SCHEMA = pl.Schema(
    {
        "qseqid": pl.String(),
        "sseqid_1": pl.String(),
        "pident_1": pl.Float64(),
        "length_1": pl.UInt32(),
        "mismatch_1": pl.UInt32(),
        "gapopen_1": pl.UInt32(),
        "qstart_1": pl.UInt32(),
        "qend_1": pl.UInt32(),
        "sstart_1": pl.UInt64(),
        "send_1": pl.UInt64(),
        "qlen_1": pl.UInt32(),
        "slen_1": pl.UInt64(),
        "evalue_1": pl.Float64(),
        "bitscore_1": pl.Float64(),
        "sseqid_2": pl.String(),
        "pident_2": pl.Float64(),
        "length_2": pl.UInt32(),
        "mismatch_2": pl.UInt32(),
        "gapopen_2": pl.UInt32(),
        "qstart_2": pl.UInt32(),
        "qend_2": pl.UInt32(),
        "sstart_2": pl.UInt64(),
        "send_2": pl.UInt64(),
        "qlen_2": pl.UInt32(),
        "slen_2": pl.UInt64(),
        "evalue_2": pl.Float64(),
        "bitscore_2": pl.Float64(),
        "total_coverage": pl.Int32(),
        "overlap_length": pl.Int32(),
        "coverage_1_only": pl.Int32(),
        "coverage_2_only": pl.Int32(),
        "chimeric": pl.Boolean(),
        "coordinate_in_1": pl.Int32(),
        "coordinate_in_2": pl.Int32(),
        "family": pl.String(),
        "species.taxid": pl.UInt64(),
        "srr": pl.String(),
    }
)

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
    """
    From a polars dataframe containing the columns:
        - slen: length of the target sequence
        - sstart: start position of the hit
        - send: end position of the hit
    Returns a polars dataframe like:
        ┌─────┬─────┬─────┬─────┬───┬──────┬──────┬──────┬──────┐
        │ 0   ┆ 1   ┆ 2   ┆ 3   ┆ … ┆ 1331 ┆ 1332 ┆ 1333 ┆ 1334 │
        │ --- ┆ --- ┆ --- ┆ --- ┆   ┆ ---  ┆ ---  ┆ ---  ┆ ---  │
        │ i64 ┆ i64 ┆ i64 ┆ i64 ┆   ┆ i64  ┆ i64  ┆ i64  ┆ i64  │
        ╞═════╪═════╪═════╪═════╪═══╪══════╪══════╪══════╪══════╡
        │ 0   ┆ 0   ┆ 0   ┆ 1   ┆ … ┆ 0    ┆ 0    ┆ 0    ┆ 0    │
        │ 0   ┆ 0   ┆ 0   ┆ 1   ┆ … ┆ 0    ┆ 0    ┆ 0    ┆ 0    │
        │ 0   ┆ 0   ┆ 0   ┆ 1   ┆ … ┆ 1    ┆ 0    ┆ 0    ┆ 0    │
        │ 0   ┆ 0   ┆ 0   ┆ 1   ┆ … ┆ 1    ┆ 0    ┆ 0    ┆ 0    │
        │ 0   ┆ 0   ┆ 0   ┆ 0   ┆ … ┆ 0    ┆ 0    ┆ 1    ┆ 1    │
        │ 0   ┆ 0   ┆ 0   ┆ 0   ┆ … ┆ 0    ┆ 0    ┆ 0    ┆ 0    │
        └─────┴─────┴─────┴─────┴───┴──────┴──────┴──────┴──────┘
    Each row represents a read, with columns indicating the position along the target sequence.
    0 and 1 represent whether the read covers the position or not.
    """
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
    """
    If start position is greater than end position, swap them.
    """
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
    """
    From a polars dataframe containing the columns:
        - slen: length of the target sequence
        - sstart: start position of the hit
        - send: end position of the hit
    Return a dataframe containing one column with coverage for each position on the target sequence, like:
        ┌──────────┐
        │ column_0 │
        │ ---      │
        │ i64      │
        ╞══════════╡
        │ 0        │
        │ 1        │
        │ 1        │
        │ 2        │
        │ 2        │
        └──────────┘
    """
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
    """
    Returns a tuple containing
    * a pandas DataFrame with coverage information
    * a boolean indicating if subsampling was needed
    The returned polars dataframe has one single column:
        ┌──────────┐
        │ column_0 │
        │ ---      │
        │ i64      │
        ╞══════════╡
        │ 0        │
        │ 1        │
        │ 1        │
        │ 2        │
        │ 2        │
        │ …        │
        │ 0        │
        │ 0        │
        │ 3        │
        │ 4        │
        │ 4        │
        └──────────┘
    It represents the coverage at each position of the target sequence.
    """

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
        chunk = lf.slice(i, CHUNKSIZE).select(["sstart", "send", "slen"]).collect()
        if chunk.is_empty():
            break
        coverage_dfs.append(get_chunk_coverage(chunk))

    return pl.concat(coverage_dfs, how="vertical"), is_subsampled


def get_series_to_plot(df: pl.DataFrame, log: bool, new_name: str) -> pd.DataFrame:
    """
    Returns a pandas DataFrame containing a single column
    If log is True, computes 10 * log10(x + 1)
    Else, returns the original column
    """
    col = df.columns[0]
    # computes 10 * log10(x + 1) if needed
    expr = 10 * (pl.col(col) + pl.lit(1)).log10() if log else pl.col(col)
    expr = expr.alias(new_name)
    return df.select(expr).to_pandas()


def get_outfile(target: str, srr_id: str, log: bool, is_subsampled: bool) -> str:
    """
    Returns the output file name for the given target, SRR ID, log status, and subsampling status.
    """
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
    log: bool,
):
    fig, ax = plt.subplots(figsize=FIGZISE)

    xlabel = f"SRR: {srr_id} / TE: {target} "

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PLOTTING NON CHIMERAS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if not non_chimera_coverage_df.is_empty():
        df = get_series_to_plot(non_chimera_coverage_df, log, "Non chimeras")
        df.plot(
            color=NON_CHIMERA_COLOR,
            ax=ax,
            xlabel=xlabel,
            **COMMON_PLOT_PARAMS,
        )
        del df

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PLOTTING CHIMERAS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # chimera_coverage is necessarily not empty because we filtered it beforehand
    df = get_series_to_plot(chimera_coverage_df, log, "Chimeras")
    df.plot(
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

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PARSE BLAST HITS
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    hit_lf = pl.scan_csv(
        args.blast_hit_file,
        separator="\t",
        schema=BLAST_OUTPUT_COL_SCHEMA,
        has_header=False,
    ).with_row_index()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # PARSE CHIMERA DATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    chimera_lf = pl.scan_csv(
        args.chimera_file,
        has_header=True,
        schema=CHIMERA_COL_SCHEMA,
    )

    if chimera_lf.limit(1).collect().is_empty():
        logger.warning("Chimera file is empty")
        sys.exit(0)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # GETTING LIST OF TARGET SEQUENCES
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # looping through target sequences and making coverage plot
    logger.info("Extracting target sequences")
    target_sequences = hit_lf.select("sseqid").unique().collect().to_series().to_list()

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MERGING BLAST HITS WITH CHIMERA DATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    logger.info("Merging blast hits with chimera data")

    # getting blast hits corresponding to chimera rows
    # choosing either one of the two pairs of coordinates:
    # for each qseqid in chimera_lf, we want to find the corresponding blast results with
    # EITHER sstart = sstart_1 & send = send_1
    # OR     sstart = sstart_2 & send = send_2
    chimera_hit_lf = (
        hit_lf.join(chimera_lf, how="full", on="qseqid", suffix="_chimera")
        .filter(
            (pl.col("qseqid") == pl.col("qseqid_chimera"))
            & (
                (
                    (pl.col("sseqid") == pl.col("sseqid_1"))
                    & (pl.col("qstart") == pl.col("qstart_1"))
                    & (pl.col("qend") == pl.col("qend_1"))
                    & (pl.col("sstart") == pl.col("sstart_1"))
                    & (pl.col("send") == pl.col("send_1"))
                )
                | (
                    (pl.col("sseqid") == pl.col("sseqid_2"))
                    & (pl.col("qstart") == pl.col("qstart_2"))
                    & (pl.col("qend") == pl.col("qend_2"))
                    & (pl.col("sstart") == pl.col("sstart_2"))
                    & (pl.col("send") == pl.col("send_2"))
                )
            )
        )
        .select(["sseqid", "sstart", "send", "slen", "index"])
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # SEPARATING BLAST HITS INTO CHIMERA AND NON-CHIMERA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    logger.info("Dividing hit table into chimera and non-chimera hits")

    chimera_indexes = (
        chimera_hit_lf.select("index").unique().collect().to_series().to_list()
    )
    non_chimera_hit_lf = hit_lf.filter(~pl.col("index").is_in(chimera_indexes))

    chimera_hit_lf = chimera_hit_lf.drop("index")
    non_chimera_hit_lf = non_chimera_hit_lf.drop("index")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MAKING COVERAGE PLOTS FOR EACH TARGET SEQUENCE
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    logger.info("Computing coverages for each target sequence and making plots")
    for target in tqdm(target_sequences):
        # getting coverage for chimera reads
        chimera_coverage_df, is_subsampled = get_coverage(
            chimera_hit_lf, target, is_chimera=True
        )

        # we don't want to plot figures for targets for which there is no chimera
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
            log=False,
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
