#!/usr/bin/env Rscript
options(error = traceback)
######################################################################################################
################## Finding chimeras in short high - throughtput sequencing reads #######################
######################################################################################################

# This script aims to identify chimeric short reads resulting from some form of recombination between two genomes. The script only searches for evidence of recombination within reads.
# It is not designed to find recombination junction between reads (in the case of paired - reads). But paired - reads can be used to search for recombination events within reads.
# The reads are first used as queries to perform a blastn search on all genomes between which recombination events are searched.
# blastn (v2.4.0+) command used to generate blastn outputs (default option megablast) for genome1 : blastn  - query reads.fasta  - db genome1.fasta  - outfmt 6  - max_target_seqs 2  - out blastOnGenome1.txt
# The same command should be used to generate blastn outputs for other genomes of interest.
# The same command should be used to generate blastn outputs for other genomes of interest.

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("arrow"))
suppressPackageStartupMessages(library("optparse"))
library(data.table)
library(optparse)
library(dplyr)
library(arrow)
library(stringr)

TMP_FOLDER <- "tmp"

#####################################################
#####################################################
# CHIMERA DETECTION FUNCTIONS
#####################################################
#####################################################

remove_overlapping_alignments <- function(blast, min_overlap_for_dropping) {
    #' among HSPs covering the same ≥20bp region of the same read, selects the one of max score

    dt = copy(blast)
    #sorts HSPs by coordinates
    setorder(dt, qstart, qend)

    if (nrow(dt) == 1) {
        return(dt)
    }

    #row(s) of HSP to discard
    to_discard = 1
    while(length(to_discard) > 0) {

        #HSPs at rows with odd numbers
        odd = 1:(nrow(dt) - 1)
        #HSPs at rows with even numbers
        even = odd + 1

        # if is TRUE for consecutive HSPs that overlap by at least 20bp
        f = abs(dt$qstart[even] - dt$qend[odd]) + 1 > min_overlap_for_dropping

        #puts consecutive overlapping HSPs in two columns of a data table (odd rows at the left)
        candidates = data.table(odds = odd[f], evens = even[f])
        #determine the one to discard according to bitscore
        to_discard = candidates[, ifelse(dt$bitscore[odds] < dt$bitscore[evens], odds, evens)]

        # checking if to_discard is not full of NA (which would result in an infinite loop)
        if ( all(is.na(to_discard)) ) {
            break
        }
        if(length(to_discard) > 0) {
            #and removes them from the blast results
            dt = dt[ -unique(to_discard),]
        }
    }

    return(dt)
}


get_best_hit_per_read <- function(dt) {
  # get the best hits based on:
  # 1: bitscore
  # 2: pident
  # 3: length
  # these three filters are used only in the case where multiple hits have the same bitscore
  best_hit_blast <- dt |>
      slice_max(bitscore, with_ties = TRUE) |> # get all hits with the highest bitscore
      slice_max(pident, with_ties = TRUE) |> # get all hits with the highest pident
      slice_max(length, with_ties = FALSE) # get only one hit (the first in the group) with the highest alignment length

  return(best_hit_blast)
}


cross_dataframes <- function(df1, df2) {
    # Make a cross product of both dataframes

    merged <- df1 |>
      mutate(.key = 1) |> # add a dummy key column to merge on
      full_join(
        df2 |> mutate(.key = 1),
        by = ".key",
        suffix = c("_1", "_2")
      ) |>
      select(-c(qseqid_2, .key)) |> # remove unnecessary columns
      rename(qseqid = qseqid_1)

    different_qlen <- merged |> filter(qlen_1 != qlen_2)
    # checking (just in case)
    if ( nrow(different_qlen) > 0 ) {
        warning("Some query read lengths do not correspond!")
    }

    return(merged)
}

get_total_coverage <- function(blast){
    # the length of the aligned region of a read (1 + 2)
    return( with(blast, qend_2 - qstart_1 + 1) )
}


get_overlap_length <- function(blast){
    #the overlap between aligned regions (on 1 and 2)
    return( with(blast, qend_1 - qstart_2 + 1) )
}

get_coverage_1_only <- function(blast){
    #part of read that only aligns to 1
    return( with(blast, qstart_2 - qstart_1) )
}

get_coverage_2_only <- function(blast){
    #part of read that only aligns on 2
    return( with(blast, qend_2 - qend_1) )
}


get_chimeric_reads <- function(dt, min_total_coverage, min_coverage_original_sequence, min_overlap, max_overlap) {
    #' finds Chimeric reads = those that partly blast on the virus and the host.
    #' The read has to partly align on the virus genome at the beginning and on the host at the end (or vice versa)
    #' with some (low) overlap between the 2 regions that align.
    #' They also must have a region where they only align on the virus and another that only align on the host.
    #' Or else, it could just reflect some contamination with host DNA that partly resemble the virus DNA

    # annotating reads as chimeric or non chimerix
    dt$chimeric <- dt$total_coverage >= min_total_coverage &
                dt$overlap_length >= min_overlap &
                dt$overlap_length <= max_overlap &
                dt$coverage_1_only >= min_coverage_original_sequence &
                dt$coverage_2_only >= min_coverage_original_sequence

    # keeping only chimeric reads
    dt = dt |> filter(chimeric == TRUE)

    return(dt)
}

compute_coverages_and_overlap <- function(dt, min_total_coverage) {
    #' Determines the overlap = length of the microhomology between parent sequences in a chimeric read
    #' using a data table (dt) of chimeric reads (output of join_best_hits_on_1_and_2_for_each_read)

    #minimum alignment length of a read against 1 and 2
    if(min_total_coverage <=  1) {
        min_total_coverage = min_total_coverage * dt$qlen_1
    } else {
        min_total_coverage = rep(min_total_coverage, nrow(dt))
    }

    total_coverage <- get_total_coverage(dt)
    overlap_length <- get_overlap_length(dt)
    coverage_1_only <- get_coverage_1_only(dt)
    coverage_2_only <- get_coverage_2_only(dt)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # corrections
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #if the overlap between the aligned parts of the read is greater than the length of aligned region, we need to swap terms
    f = abs(overlap_length) > abs(total_coverage)
    f[is.na(f)] = FALSE

    #temporary vector for the swap
    tmp = overlap_length
    overlap_length[f] = total_coverage[f]
    total_coverage[f] = tmp[f]

    tmp = coverage_1_only
    coverage_1_only[f] = -coverage_2_only[f]
    coverage_2_only[f] = -tmp[f]

    #could be negative if the alignement on 1 is totally included in that on 2
    coverage_2_only[coverage_2_only < 0] = 0
    coverage_1_only[coverage_1_only < 0] = 0
    overlap_length[f] = total_coverage[f]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    dt$total_coverage <- total_coverage
    dt$overlap_length <- overlap_length
    dt$coverage_1_only <- coverage_1_only
    dt$coverage_2_only <- coverage_2_only

    return(dt)
}


get_coordinate_in_1 <- function(dt, dec) {
  #' finds coordinates of host sequences in 1
  coordH = rep(NA, nrow(dt))

  f = with(dt, chimeric & qstart_2 < qstart_1 & sstart_2 < send_2) #alignment with host at the beginning of read, "+" direction
  coordH[f] = dt$send_2[f]

  f = with(dt, chimeric & qstart_2 > qstart_1 & sstart_2 > send_2) #alignment at the end of read, direction minus ( - )
  coordH[f]  = dt$sstart_2[f]  #same as above

  f = with(dt, chimeric & qstart_2 < qstart_1 & sstart_2 > send_2)  #alignment at beginning of read, minus direction
  coordH[f]  = dt$send_2[f] - 1	#here the insertion site is the start of match in the sseqid, minus 1 for consistency (as we use the base BEFORE the insertion site)

  f = with(dt, chimeric & qstart_2 > qstart_1 & sstart_2 < send_2) #alignment at end of read, + direction
  coordH[f]  = dt$sstart_2[f] - 1   #same as above

  return(coordH)
}


get_coordinate_in_2 <- function(dt) {
  #' finds insertion coordinates of host sequences in the virus genome.
  coord = rep(NA, nrow(dt))

  f = with(dt, chimeric & qstart_1 < qstart_2 & sstart_1 < send_1) #alignment with virus at the beginning of read, "+" direction
  coord[f] = dt$send_1[f]  #the insertion site is the end of match in the virus (hence largest coordinate in sseqid)

  f2 = with(dt, chimeric & qstart_1 > qstart_2 & sstart_1 > send_1) #alignment at the end of read (101), direction minus ( - )
  coord[f2] = dt$sstart_1[f2]  #same as above

  f3 = with(dt, chimeric & qstart_1 < qstart_2 & sstart_1 > send_1)  #match at beginning of read, minus direction
  coord[f3] = dt$send_1[f3] - 1  #here the insertion site is the start of match in the sseqid, minus 1 for consistency (as we use the base BEFORE the insertion site)

  f4 = with(dt, chimeric & qstart_1 > qstart_2 & sstart_1 < send_1) #match at end of read, + direction
  coord[f4] = dt$sstart_1[f4] - 1   #same as above

  return(coord)
}


find_chimeras <- function (dt1, dt2) {

    # CONSTANTS
    MIN_OVERLAP_FOR_DROPPING <- 20
    MIN_TOTAL_COVERAGE <- 0.9
    MIN_OVERLAP <- -5
    MAX_OVERLAP <- 20
    MIN_COVERAGE_ORIGINAL_SEQUENCE <- 16

    # remove alignments within the blast object that are overlapping for a same read, keeping the best alignment (best bitscore)
    dt1 <- remove_overlapping_alignments(dt1, MIN_OVERLAP_FOR_DROPPING)
    dt2 <- remove_overlapping_alignments(dt2, MIN_OVERLAP_FOR_DROPPING)

    dt1 <- get_best_hit_per_read(dt1)
    dt2 <- get_best_hit_per_read(dt2)

    dt <- cross_dataframes(dt1, dt2)

    dt <- compute_coverages_and_overlap(dt, MIN_TOTAL_COVERAGE)

    # find chimeric reads. Four parameters can be set:
    # 1 - proportion of the read which is aligned, cumulating alignment length on the 2 genomes (here 0.9, meaning 90% of the read has to be aligned)
    # 2 - maximum number of bases inserted between the 2 genomes at recombination point, reflecting non - templated nucleotide additions (here 5)
    # 3 - maximum overlap in the alignment with the 2 genomes at the recombination points, reflects the presence of homology between the 2 genomes at the recombination point (here 20)
    # 4 - minimum alignment length on one genome only. Here the read has to be aligned over at least 16 bp on the genome 1 only and over at least 16 bp on genome 2 only
    dt = get_chimeric_reads(
        dt,
        MIN_TOTAL_COVERAGE,
        MIN_COVERAGE_ORIGINAL_SEQUENCE,
        MIN_OVERLAP,
        MAX_OVERLAP
    )

    # insert overlap column containing the number of nucleotides shared between the 2 genomes at the recombination point
    dt$overlap_length = get_overlap_length(dt)

    # insert column containing coordinate of the recombination point in 1
    dt$coordinate_in_1 = get_coordinate_in_1(dt)

    # insert column containing coordinate of the recombination point in 2
    dt$coordinate_in_2 = get_coordinate_in_2(dt)

    return(dt)

}

#####################################################
#####################################################
# PARSING AND DATA MANAGEMENT FUNCTIONS
#####################################################
#####################################################


get_args <- function() {

    option_list <- list(
        make_option("--hits-1", dest = 'blast_hits_1_file', help = "Path to first Blast output file"),
        make_option("--hits-2", dest = 'blast_hits_2_file', help = "Path to second Blast output file for genome"),
        make_option("--family", help = "Family name / Taxid"),
        make_option("--species", help = "Species Taxid"),
        make_option("--srr", dest = 'srr_id', help = "ID of the corresponding SRR"),
        make_option("--out", dest = 'outfile', help = "Output file name")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Find chimeric sequences"
        ))

    return(args)
}

get_blast_output_schema <- function() {
  # Sets the schema for the blast output dataframe
  # qseqid      query or source (gene) sequence id
  # sseqid      sseqid or target (reference genome) sequence id
  # pident      percentage of identical positions
  # length      alignment length (sequence overlap)
  # mismatch    number of mismatches
  # gapopen     number of gap openings
  # qstart      start of alignment in query
  # qend        end of alignment in query
  # sstart      start of alignment in sseqid
  # send        end of alignment in sseqid
  # qlen        [additional column] query sequence length
  # slen        [additional column] sseqid sequence length
  # evalue      expect value
  # bitscore    bitscore
  #
  # see https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
  # columns defined in conf/modules/blast.config
  arrow::schema(
    qseqid     = utf8(),
    sseqid     = utf8(),
    pident     = float64(),
    length     = int32(),
    mismatch   = int32(),
    gapopen    = int32(),
    qstart     = int32(),
    qend       = int32(),
    sstart     = int32(),
    send       = int32(),
    qlen       = int32(),
    slen       = int32(),
    evalue     = float64(),
    bitscore   = float64()
  )
}

get_partitioned_dataset_schema <- function() {
  new_schema <- arrow::schema(bucket = utf8())
  arrow::schema(c(get_blast_output_schema()$fields, new_schema$fields))
}

get_arrow_dataset_blast_output <- function(file) {
  arrow::open_dataset(
    file,
    schema = get_blast_output_schema(),
    format = "text",
    delimiter = "\t"
  )
}

get_partitioned_dataset <- function(file) {
  arrow::open_dataset(file, schema = get_partitioned_dataset_schema())
}

get_bucket_id <- function(col) {
  stringr::str_sub(as.character(col), -3, -1)
}

partition_dataset_on_qseqid <- function(ds, outfile) {
  ds |>
    mutate(bucket = as.character(get_bucket_id(qseqid))) |> # attributes an integer between 0 and 999 to each qseqid based on last 3 digits
    write_dataset(outfile, partitioning = "bucket") # writes one partition for each bucket
}

get_unique_qseqids <- function(ds) {
  ds |> distinct(qseqid) |> pull(qseqid, as_vector = TRUE)
}

get_common_qseqids <- function(ds1, ds2) {
  uniques_ds1_qseqids <- get_unique_qseqids(ds1)
  uniques_ds2_qseqids <- get_unique_qseqids(ds2)
  intersect(uniques_ds1_qseqids, uniques_ds2_qseqids)
}

get_dataframe_subset <- function(ds, bucket_id, read_id) {
  ds |> filter(bucket == bucket_id, qseqid == read_id) |> select(-bucket) |> collect()
}

process_batches <- function(ds1, ds2) {

  common_qseqids <- get_common_qseqids(ds1, ds2)
  nb_unique_qseqids <- length(common_qseqids)

  nb_processed <- 0
  for (read_id in common_qseqids) {

    bucket_id <- get_bucket_id(read_id)
    df1 <- get_dataframe_subset(ds1, bucket_id, read_id)
    df2 <- get_dataframe_subset(ds2, bucket_id, read_id)

    chimera_dt <- find_chimeras(df1, df2)

    tmp_outfile <- paste0(TMP_FOLDER, "/", read_id, ".csv")
    write.table(chimera_dt, tmp_outfile, sep = ',', row.names = FALSE, quote = FALSE)

    nb_processed <- nb_processed + 1
    pct_done <- 100 * nb_processed / nb_unique_qseqids
    message(paste0("Qseqid ", read_id, " done. ", nb_processed, " qseqid processed (",  format(round(pct_done, 2), nsmall = 2), "% of total)."))
  }

  message("Collecting all tempopary chimera results into one single dataframe")
  chimera_filenames <- dir(TMP_FOLDER, pattern = "\\.csv$")

  if (length(chimera_filenames) == 0) {
    message("No CSV files found.")
    df <- data.frame()
  } else {
    csv_files <- paste0(TMP_FOLDER, "/", chimera_filenames)
    df <- do.call(rbind, lapply(csv_files, read.csv))
  }

  return(df)
}

add_metadata <- function(dt, family, species, srr) {
    message("Adding metadata")
    dt$family = family
    dt$species.taxid = species
    dt$srr = srr
    return(dt)
}

export_data <- function(dt, filename) {
    message(paste('Exporting data to:', filename, "\n"))
    write.table(dt, filename, sep = ',', row.names = FALSE, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

args <- get_args()

if ( file.size(args$blast_hits_1_file) > 0 && file.size(args$blast_hits_2_file) > 0 ) {

    # parse and write with partitioning on qesqid, so that subsequent filtering get much faster
    message("Parsing dataset 1")
    ds1 <- get_arrow_dataset_blast_output(args$blast_hits_1_file)
    message("Partitioning dataset 1...")
    partition_dataset_on_qseqid(ds1, "ds1_partitioned")
    rm(ds1)

    message("Parsing dataset 2")
    ds2 <- get_arrow_dataset_blast_output(args$blast_hits_2_file)
    message("Partitioning dataset 2...")
    partition_dataset_on_qseqid(ds2, "ds2_partitioned")
    rm(ds2)

    # read the partitioned datasets
    ds1 <- get_partitioned_dataset("ds1_partitioned")
    ds2 <- get_partitioned_dataset("ds2_partitioned")

    ds1_nrows <- nrow(ds1)
    ds2_nrows <- nrow(ds2)
    message(paste("File 1 dataset has", ds1_nrows, "rows."))
    message(paste("File 2 dataset has", ds2_nrows, "rows."))

    dir.create(TMP_FOLDER)
    df <- process_batches(ds1, ds2)

    message("Removing partitioned datasets...")
    unlink("ds1_partitioned", recursive = TRUE)
    unlink("ds2_partitioned", recursive = TRUE)

} else {
    message("At least one input file is empty")
    df <- data.table()
}

# we want to export files anyway, even when there is no read
# this helps to keep track of SRRs that have been processed until the end
nb_chimeras <- nrow(df)
if ( nb_chimeras == 0 ) {
    message("\nNo chimeras found")
    file.create(args$outfile)
} else {
    message(paste("\nFound ", nb_chimeras, " chimeras in total."))
    df <- add_metadata(df, args$family, args$species, args$srr)
    export_data(df, args$outfile)
}
