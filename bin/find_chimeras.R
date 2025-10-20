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
library(data.table)
library(optparse)
library(dplyr)

#####################################################
#####################################################
# FUNCTIONS
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

# functions used for the search and analysis of chimeric reads

remove_overlapping_alignments <- function(blast, min_overlap_for_dropping) {
    #' among HSPs covering the same ≥20bp region of the same read, selects the one of max score

    dt = copy(blast)
    #sorts HSPs by coordinates
    setorder(dt, qseqid, qstart, qend)

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

        #f is TRUE for consecutive HSPs that overlap by at least 20bp
        f = abs(dt$qstart[even] - dt$qend[odd]) + 1 > min_overlap_for_dropping &
            dt$qseqid[odd] ==  dt$qseqid[even]

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
  # 3: qlen
  # these three filters are used only in the case where multiple hits have the same bitscore
  best_hit_blast <- dt %>%
      group_by(qseqid) %>%
      slice_max(bitscore, with_ties = TRUE) %>% # get all hits with the highest bitscore
      slice_max(pident, with_ties = TRUE) %>% # get all hits with the highest pident
      slice_max(length, with_ties = FALSE) %>% # get only one hit (the first in the group) with the highest alignment length
      ungroup()

  return(best_hit_blast)
}


get_reads_with_hits_on_both <- function(dt1, dt2) {
    # Join both dataframes (inner join) to keep only reads with hits on both target and genome

    #merges these tables (all = FALSE: inner join)
    merged <- merge(dt1, dt2, by = "qseqid", all = FALSE, suffixes = c("_1","_2"))
    message(paste("Obtained", nrow(merged), "reads having hits on both 1 and 2"))

    # checking (just in case)
    if ( any(merged$qlen_1 != merged$qlen_2) ) {
        error("Read lengths do not correspond")
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
    dt = dt %>% filter(chimeric == TRUE)

    return(dt)
}

compute_coverages_and_overlap <- function(dt, min_total_coverage) {
    #' Determines the overlap = length of the microhomology between parent sequences in a chimeric read
    #' using a data table (dt) of chimeric reads (output of join_best_hits_on_1_and_2_for_each_read)

    #minimum alignment length of a read against 1 and 2
    if(min_total_coverage <=  1) {
        min_total_coverage = min_total_coverage * dt$qlen
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
    message("Removing overlapping alignments for each read")
    dt1 <- remove_overlapping_alignments(dt1, MIN_OVERLAP_FOR_DROPPING)
    dt2 <- remove_overlapping_alignments(dt2, MIN_OVERLAP_FOR_DROPPING)

    message("Keeping best hit per read")
    dt1 <- get_best_hit_per_read(dt1)
    dt2 <- get_best_hit_per_read(dt2)

    message("Keeping reads having a hit on both target and genome")
    dt <- get_reads_with_hits_on_both(dt1, dt2)

    message("Computing coverages and overlap")
    dt <- compute_coverages_and_overlap(dt, MIN_TOTAL_COVERAGE)

    # find chimeric reads. Four parameters can be set:
    # 1 - proportion of the read which is aligned, cumulating alignment length on the 2 genomes (here 0.9, meaning 90% of the read has to be aligned)
    # 2 - maximum number of bases inserted between the 2 genomes at recombination point, reflecting non - templated nucleotide additions (here 5)
    # 3 - maximum overlap in the alignment with the 2 genomes at the recombination points, reflects the presence of homology between the 2 genomes at the recombination point (here 20)
    # 4 - minimum alignment length on one genome only. Here the read has to be aligned over at least 16 bp on the genome 1 only and over at least 16 bp on genome 2 only
    message("Detecting chimeric reads...")
    dt = get_chimeric_reads(
        dt,
        MIN_TOTAL_COVERAGE,
        MIN_COVERAGE_ORIGINAL_SEQUENCE,
        MIN_OVERLAP,
        MAX_OVERLAP
    )

    # insert overlap column containing the number of nucleotides shared between the 2 genomes at the recombination point
    message("Adding length of overlaps between 1 and 2")
    dt$overlap_length = get_overlap_length(dt)

    # insert column containing coordinate of the recombination point in 1
    dt$coordinate_in_1 = get_coordinate_in_1(dt)

    # insert column containing coordinate of the recombination point in 2
    dt$coordinate_in_2 = get_coordinate_in_2(dt)

    return(dt)

}

parse_blast_hit_file <- function(blast_hits_file) {
    #' Parse blast hit result file
    # Add propre header to dataframes
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

    blast_hits_dt = fread(blast_hits_file)
    if ( nrow(blast_hits_dt) == 0 ) {
        return(blast_hits_dt)
    }
    # columns defined in conf/modules/blast.config
    BLAST_OUTFMT6_COLS <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "qlen", "slen", "evalue", "bitscore")
    setnames(blast_hits_dt, BLAST_OUTFMT6_COLS)
    return(blast_hits_dt)
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

blast_hits_1_dt <- parse_blast_hit_file(args$blast_hits_1_file)
blast_hits_2_dt <- parse_blast_hit_file(args$blast_hits_2_file)

if ( nrow(blast_hits_1_dt) == 0 || nrow(blast_hits_2_dt) == 0 ) {
    warning("At least one input file is empty")
    chimera_dt <- data.table()
} else {
    chimera_dt <- find_chimeras(blast_hits_1_dt, blast_hits_2_dt)
}

# we want to export files anyway, even when there is no read
# this helps t keep tracjk of SRRs that have been processed until the end
nb_chimeras <- nrow(chimera_dt)
if ( nb_chimeras == 0 ) {
    message("\nNo chimeras found")
    file.create(args$outfile)
    file.create("chimeric_reads.txt")
} else {
    message(paste("\nFound ", nb_chimeras, " chimeras"))
    # writing the names of the chimeric reads in a file
    write(chimera_dt$qseqid, file = "chimeric_reads.txt", ncolumns = 1)
    chimera_dt <- add_metadata(chimera_dt, args$family, args$species, args$srr)
    export_data(chimera_dt, args$outfile)
}
