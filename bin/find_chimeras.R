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
library(data.table)
library(optparse)

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

remove_overlapping_alignments = function(blast, min_overlap_for_dropping) {
    #' among HSPs covering the same ≥20bp region of the same read, selects the one of higher evalue
    
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
        #determine the one to discard according to evalue	
        to_discard = candidates[, ifelse(dt$evalue[odds] < dt$evalue[evens], odds, evens)]					
        
        # checking if to_discard is full of NA (resulting in infinite loop)
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


join_best_hits_on_1_and_2_for_each_read = function(blast) {
    #' Puts the 2 best HSPs on the same read (which are in different rows) at the same row in a new table,
    #' so that it can be determined  afterwards if such reads is chimeric

    #randomise the order of HSPs in the blast results table
    blast = blast[sample(1:nrow(blast))]
    #selects reads that have several HSPs
    blast = blast[duplicated(qseqid) | duplicated(qseqid, fromLast = T)]
    #puts the best HSPs on top
    setorder(blast, qseqid, -evalue)

    #puts best HSP on each read in a new table
    blast1 = blast[!duplicated(qseqid)]
    #puts 2nd best HSP on each read in another table
    blast2 = blast[duplicated(qseqid)]
    blast2 = blast2[!duplicated(qseqid)]

    #merges these tables (full outer join)
    merged <- merge(blast1, blast2, by = "qseqid", all = T, suffixes = c("_1","_2"))
    message(paste("Obtained", nrow(merged), "reads having hits on both 1 and 2"))

    # checking
    if ( any(merged$qlen_1 != merged$qlen_2) ) {
        error("Read lengths do not correspond")
    }

    # keep only rows where subject types are different
    merged <- merged[merged$subject_type_1 != merged$subject_type_2]
    message(paste("Kept", nrow(merged), "rows corresponding to hits on both 1 and 2"))

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


get_chimeric_state <- function(df, min_total_coverage, min_coverage_original_sequence, min_overlap, max_overlap) {
    #' finds Chimeric reads = those that partly blast on the virus and the host.
    #' The read has to partly align on the virus genome at the beginning and on the host at the end (or vice versa)
    #' with some (low) overlap between the 2 regions that align.
    #' They also must have a region where they only align on the virus and another that only align on the host.
    #' Or else, it could just reflect some contamination with host DNA that partly resemble the virus DNA

    chimeric <- df$total_coverage >= min_total_coverage &
                df$overlap_length >= min_overlap &
                df$overlap_length <= max_overlap &
                df$coverage_1_only >= min_coverage_original_sequence &
                df$coverage_2_only >= min_coverage_original_sequence

    # replacing NA by False
    chimeric[is.na(chimeric)] = F

    return(chimeric)
}

compute_coverages_and_overlap <- function(df, min_total_coverage) {
    #' Determines the overlap = length of the microhomology between parent sequences in a chimeric read
    #' using a data table (dt) of chimeric reads (output of join_best_hits_on_1_and_2_for_each_read)

    #minimum alignment length of a read against 1 and 2
    if(min_total_coverage <=  1) {
        min_total_coverage = min_total_coverage * df$qlen
    } else {
        min_total_coverage = rep(min_total_coverage, nrow(df))
    }

    total_coverage <- get_total_coverage(df)
    overlap_length <- get_overlap_length(df)
    coverage_1_only <- get_coverage_1_only(df)
    coverage_2_only <- get_coverage_2_only(df)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # corrections
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #if the overlap between the aligned parts of the read is greater than the length of aligned region, we need to swap terms
    f = abs(overlap_length) > abs(total_coverage)
    f[is.na(f)] = F

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

    df$total_coverage <- total_coverage
    df$overlap_length <- overlap_length
    df$coverage_1_only <- coverage_1_only
    df$coverage_2_only <- coverage_2_only

    return(df)
}


get_coordinate_in_1 = function(dt, dec) {
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


get_coordinate_in_2 = function(dt) {
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


find_chimeras <- function (blast_hits_1_df, blast_hits_2_df) {

    # CONSTANTS
    MIN_OVERLAP_FOR_DROPPING <- 20
    MIN_TOTAL_COVERAGE <- 0.9
    MIN_OVERLAP <- -5
    MAX_OVERLAP <- 20
    MIN_COVERAGE_ORIGINAL_SEQUENCE <- 16

    # remove alignments within the blast object that are overlapping for a same read, keeping the best evalue alignment
    message("Removing overlapping alignments in Blast hits against 1...")
    blast_hits_1_df = remove_overlapping_alignments(blast_hits_1_df, MIN_OVERLAP_FOR_DROPPING)
    message("Removing overlapping alignments in blast hits against 2...")
    blast_hits_2_df = remove_overlapping_alignments(blast_hits_2_df, MIN_OVERLAP_FOR_DROPPING)

    message("Concatenating both Blast hit dataframes...")
    hits_df = rbind(blast_hits_1_df, blast_hits_2_df)

    message("Removing overlapping alignments in merged Blast hits...")
    hits_df = remove_overlapping_alignments(hits_df, MIN_OVERLAP_FOR_DROPPING)

    # blastn_noOverlap is sorted according to 'qseqid' and 'evalue' columns
    setorder(hits_df, qseqid, -evalue)

    # for a given read only the 2 best - evalue hits are kept. The description lines of these two hits are then merged on the same row
    message("Keeping the 2 best hits for each read...")
    hits_df = join_best_hits_on_1_and_2_for_each_read(hits_df)

    message("Computing coverages and overlap")
    hits_df <- compute_coverages_and_overlap(hits_df, MIN_TOTAL_COVERAGE)

    # find chimeric reads. Four parameters can be set:
    # 1 - proportion of the read which is aligned, cumulating alignment length on the 2 genomes (here 0.9, meaning 90% of the read has to be aligned)
    # 2 - maximum number of bases inserted between the 2 genomes at recombination point, reflecting non - templated nucleotide additions (here 5)
    # 3 - maximum overlap in the alignment with the 2 genomes at the recombination points, reflects the presence of homology between the 2 genomes at the recombination point (here 20)
    # 4 - minimum alignment length on one genome only. Here the read has to be aligned over at least 16 bp on the genome 1 only and over at least 16 bp on genome 2 only
    message("Detecting chimeric reads...")
    hits_df$chimeric = get_chimeric_state(
        hits_df,
        MIN_TOTAL_COVERAGE,
        MIN_COVERAGE_ORIGINAL_SEQUENCE,
        MIN_OVERLAP,
        MAX_OVERLAP
    )
    
    # generate data table containing only chimeric reads
    message("Keeping only chimeric reads")
    chimera_df = hits_df[chimeric == T]

    # insert overlap column containing the number of nucleotides shared between the 2 genomes at the recombination point
    message("Adding length of overlaps between 1 and 2")
    chimera_df$overlap_length = get_overlap_length(chimera_df)

    # insert column showing which chimeric reads may be PCR duplicates, i.e. which chimeric reads have identical alignement coordinates on both genomes
    chimera_df[, pcr_duplicate:= paste(qstart_1 - sstart_1, qend_1 - send_1, sseqid_1, qstart_2 - sstart_2, qend_2 - send_2, sseqid_2)]

    # insert column containing coordinate of the recombination point in 1
    chimera_df$coordinate_in_1 = get_coordinate_in_1(chimera_df)

    # insert column containing coordinate of the recombination point in 2
    chimera_df$coordinate_in_2 = get_coordinate_in_2(chimera_df)

    # insert column showing the respective orientation (same or opposite) of the sequences involved in intra - genome chimeras
    chimera_df[, inverted:= sign(send_1 - sstart_1) != sign(send_2 - sstart_2)]

    return(chimera_df)

}

parse_blast_hit_file <- function(blast_hits_file, sseqid_type) {
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
    # qlen        query sequence length
    # slen        sseqid sequence length
    # evalue      expect value
    # bitevalue    bit evalue
    #
    # see https://www.metagenomics.wiki/tools/blast/blastn-output-format-6

    blast_hits_df = fread(blast_hits_file)
    if ( nrow(blast_hits_df) == 0 ) {
        return(blast_hits_df)
    }
    # columns defined in conf/modules/blast.config
    BLAST_OUTFMT6_COLS <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "qlen", "slen", "evalue", "bitevalue")
    setnames(blast_hits_df, BLAST_OUTFMT6_COLS)
    # setting sseqid type
    blast_hits_df$subject_type = sseqid_type
    return(blast_hits_df)
}

add_metadata <- function(df, family, species, srr) {
    message("Adding metadata")
    df$family = family
    df$species.taxid = species
    df$srr = srr
    return(df)
}

export_data <- function(df, filename) {
    message(paste('Exporting data to:', filename, "\n"))
    write.table(df, filename, sep = ',', row.names = FALSE, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


args <- get_args()

blast_hits_1_df <- parse_blast_hit_file(args$blast_hits_1_file, 1)
blast_hits_2_df <- parse_blast_hit_file(args$blast_hits_2_file, 2)

if ( nrow(blast_hits_1_df) == 0 || nrow(blast_hits_2_df) == 0 ) {
    warning("At least one input file is empty")
    chimera_df <- data.table()
} else {
    chimera_df <- find_chimeras(blast_hits_1_df, blast_hits_2_df)
}

# we want to export files anyway, even when there is no read
# this helps t keep tracjk of SRRs that have been processed until the end
nb_chimeras <- nrow(chimera_df)
if ( nb_chimeras == 0 ) {
    message("\nNo chimeras found")
    file.create(args$outfile)
    file.create("chimeric_reads.txt")
} else {
    message(paste("\nFound ", nb_chimeras, " chimeras"))
    # writing the names of the chimeric reads in a file
    write(chimera_df$qseqid, file = "chimeric_reads.txt", ncolumns = 1)
    chimera_df <- add_metadata(chimera_df, args$family, args$species, args$srr)
    export_data(chimera_df, args$outfile)
}
