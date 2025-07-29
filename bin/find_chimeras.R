#!/usr/bin/env Rscript

######################################################################################################
################## Finding chimeras in short high - throughtput sequencing reads #######################
######################################################################################################

# This script aims to identify chimeric short reads resulting from some form of recombination between two genomes. The script only searches for evidence of recombination within reads.
# It is not designed to find recombination junction between reads (in the case of paired - reads). But paired - reads can be used to search for recombination events within reads.
# The reads are first used as queries to perform a blastn search on all genomes between which recombination events are searched.
# blastn (v2.4.0+) command used to generate blastn outputs (default option megablast) for genome1 : blastn  - query reads.fasta  - db genome1.fasta  - outfmt 6  - max_target_seqs 2  - out blastOnGenome1.txt
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
        make_option("--target-hits", dest = 'blast_hits_1', help = "Path to Blast output file for target database"),
        make_option("--assembly-hits", dest = 'blast_hits_2', help = "Path to Blast output file for assembly"),
        make_option("--out", dest = 'outfile', help = "Output file name")
    )

    args <- parse_args(OptionParser(
        option_list = option_list,
        description = "Find chimeric sequences"
        ))

    return(args)
}

# functions used for the search and analysis of chimeric reads

removeOverlapingAlignments = function(blast, overlap = 20) {	#among HSPs covering the same ≥20bp region of the same read, selects the one of higher score
  dt = copy(blast)
  setorder(dt, readName, qStart, qEnd)							#sorts HSPs by coordinates
  toDiscard = 1													#row(s) of HSP to discard
  while(length(toDiscard)>0) {
    odd = 1:(nrow(dt) - 1) 										#HSPs at rows with odd numbers
    even = odd + 1 												#HSPs at rows with even numbers
    f = abs(dt$qStart[even] - dt$qEnd[odd]) + 1 > overlap & dt$readName[odd] ==  dt$readName[even]		#f is TRUE for consecutive HSPs that overlap by at least 20bp
    candidates = data.table(odds = odd[f], evens = even[f])											#puts consecutive overlapping HSPs in two columns of a data table (odd rows at the left)
    toDiscard = candidates[,ifelse(dt$score[odds] < dt$score[evens], odds, evens)]					#determine the one to discard according to score
    if(length(toDiscard)>0) dt = dt[ - unique(toDiscard),]											#and removes them from the blast results
    cat("*")
  }
  dt
}


mergeOneBlast = function(blast) {						#puts the 2 best HSPs on the same read (which are in different rows) at the same row in a new table, so that it can be determined  afterwards if such reads is chimeric
  blast = blast[sample(1:nrow(blast))]					#randomise the order of HSPs in the blast results table
  blast = blast[duplicated(readName) | duplicated(readName,fromLast = T)]			#selects reads that have several HSPs
  setorder(blast, readName,  - score)						#puts the best HSPs on top
  blast1 = blast[!duplicated(readName)]					#puts best HSP on each read in a new table
  blast2 = blast[duplicated(readName)]
  blast2 = blast2[!duplicated(readName)]				#puts 2nd best HSP on each read in another table
  merge(blast1, blast2,by =  "readName", all = T,suffixes = c("",".s"))			#merges these tables
}


chimericReads = function(blast,minCov,minOverlap,maxOverlap,minOnly) { #finds Chimeric reads = those that partly blast on the virus and the host. The read has to partly align on the virus genome at the beginning and on the host at the end (or vice versa) with some  (low) overlap between the 2 regions that align. They also must have a region where they only align on the virus and another that only align on the host. Or else, it could just reflect some contamination with host DNA that partly resemble the virus DNA
  if(minCov <=  1) {
    minCov = minCov * blast$readLength
  } else {
    minCov = rep(minCov,nrow(blast))    				#minimum alignment length of a read against virus and host
  }
  cov = with(blast,qEnd.s - qStart + 1)    					#the length of the aligned region of a read (virus + host)
  overlap = with(blast,qEnd - qStart.s + 1)  				#the overlap between aligned regions (on virus and host)
  vOnly = with(blast, qStart.s - qStart)   				#part of read that only aligns to the virus
  hOnly = with(blast, qEnd.s - qEnd)       				#part of read that only aligns on the host
  f = abs(overlap) > abs(cov)   						#if the overlap between the aligned parts of the read is greater than the length of aligned region, we need to swap terms
  f[is.na(f)] = F
  t = overlap       									#temporary vector for the swap
  overlap[f] = cov[f]
  cov[f] = t[f]
  t = vOnly
  vOnly[f] =  - hOnly[f]
  hOnly[f] =  - t[f]
  hOnly[hOnly < 0] = 0   									#could be negative if the alignement on the virus is totally included in that on the host
  vOnly[vOnly < 0] = 0
  chimeric = cov >=  minCov & overlap >=  minOverlap & overlap <=  maxOverlap & vOnly >=  minOnly & hOnly >=  minOnly
  chimeric[is.na(chimeric)] = F
  return(chimeric)
}

overlap = function(dt) {								#determines the overlap = length of the microhomology between parent sequences in a chimeric read, using a data table (dt) of chimeric reads (output of mergeOneBlast())
  cov = with(dt,qEnd.s - qStart + 1)        				#the length of the aligned region of a read (virus + host)
  overlap = with(dt,qEnd - qStart.s + 1) 			 		#the overlap between aligned regions (on virus and host)
  f = abs(overlap) > abs(cov)   						#if the overlap between the aligned parts of the read is greater than the length of aligned region, we need to swap terms
  f[is.na(f)] = F
  t = overlap   										#temporary vector for the swap
  overlap[f] = cov[f]
  overlap[!dt$chimeric] = NA
  overlap
}


chimericPoint = function(dt,dec = 30) {
  coordH = rep(NA,nrow(dt))
  dt$n = 1:nrow(dt)
  if (any(dt$interChimera,na.rm = T)) {
    plus = with(dt, !is.na(subject.s) & sEnd.s > sStart.s & interChimera)
    minus = with(dt, !is.na(subject.s) & sEnd.s < sStart.s & interChimera)
    coordH[plus] = dt$sEnd.s[plus] + dec
    coordH[minus] = dt$sEnd.s[minus] - dec
    warning("note that estimates for junctions between reads are not accurate")
  }

  f = with(dt, chimeric & qStart.s < qStart & sStart.s < sEnd.s) #alignment with host at the beginning of read, "+" direction
  coordH[f] = dt$sEnd.s[f]

  f = with(dt, chimeric & qStart.s > qStart & sStart.s > sEnd.s) #alignment at the end of read, direction minus ( - )
  coordH[f]  = dt$sStart.s[f]  #same as above

  f = with(dt, chimeric & qStart.s < qStart & sStart.s > sEnd.s)  #alignment at beginning of read, minus direction
  coordH[f]  = dt$sEnd.s[f] - 1	#here the insertion site is the start of match in the subject, minus 1 for consistency (as we use the base BEFORE the insertion site)

  f = with(dt, chimeric & qStart.s > qStart & sStart.s < sEnd.s) #alignment at end of read, + direction
  coordH[f]  = dt$sStart.s[f] - 1   #same as above
  return(coordH)
}


insertionCoord = function(dt, dec = 30, withOv = F) {  #finds insertion coordinates of host sequences in the virus genome. Dec is a guesstimation of half the insert size, to roughtly locate the junction point when it occured between paired reads
  dt$n = 1:nrow(dt)
  coord = rep(NA,nrow(dt))
  if (any(dt$interChimera)) {
    plus = with(dt, !is.na(subject) & sEnd > sStart & interChimera)
    minus = with(dt, !is.na(subject) & sEnd < sStart & interChimera)
    coord[plus] = dt$sEnd[plus] + dec
    coord[minus] = dt$sEnd[minus] - dec
    #warning("note that estimates for junctions between reads are not accurate")
  }
  f = with(dt, chimeric & qStart < qStart.s & sStart < sEnd) #alignment with virus at the beginning of read, "+" direction
  coord[f] = dt$sEnd[f]  #the insertion site is the end of match in the virus (hence largest coordinate in subject)

  f2 = with(dt, chimeric & qStart > qStart.s & sStart > sEnd) #alignment at the end of read (101), direction minus ( - )
  coord[f2] = dt$sStart[f2]  #same as above

  f3 = with(dt, chimeric & qStart < qStart.s & sStart > sEnd)  #match at beginning of read, minus direction
  coord[f3] = dt$sEnd[f3] - 1  #here the insertion site is the start of match in the subject, minus 1 for consistency (as we use the base BEFORE the insertion site)

  f4 = with(dt, chimeric & qStart > qStart.s & sStart < sEnd) #match at end of read, + direction
  coord[f4] = dt$sStart[f4] - 1   #same as above

  if (withOv) {
    ov = overlap(dt)
    ov[ov < 0] = 0
    coord[f | f2] = coord[f | f2] - ov[f | f2]
    coord[f3 | f4] = coord[f3 | f4] + ov[f3 | f4]
  }
  coord
}


find_chimeras <- function (blast1, blast2) {

    cols = c("readName", "subject", "identity", "alength", "mismatch", "indel", "qStart", "qEnd", "sStart", "sEnd", "score", "readLength", "Sp")

    blast1$Sp = "target"
    setnames(blast1, cols)

    blast2$Sp = "assembly"
    setnames(blast2, cols)

    # concatenate the output blast files into one data table

    CATblastn = rbind(blast1, blast2)
    setnames(CATblastn, c("readName", "subject", "identity", "alength", "mismatch", "indel", "qStart", "qEnd", "sStart", "sEnd", "score", "readLength", "sample"))


    # remove alignments within the blast object that are overlaping for a same read, keeping the best score alignment
    blast1nooverlap = removeOverlapingAlignments(blast1)
    blast2nooverlap = removeOverlapingAlignments(blast2)

    CATblastn = rbind(blast1nooverlap, blast2nooverlap)
    CATblastn_noOverlap = removeOverlapingAlignments(CATblastn)


    # blastn_noOverlap is sorted according to 'readName' and 'score' columns
    setorder(CATblastn_noOverlap, readName,  - score)

    # for a given read only the 2 best - score hits are kept. The description lines of these two hits are then merged on the same row
    mCATblastn_noOverlap = mergeOneBlast(CATblastn_noOverlap)

    # insert column containing the length of each read. In this example all reads are 151 bp but reads can be of different lengths.
    mCATblastn_noOverlap$readLength = 150

    # find chimeric reads. Four parameters can be set:
    # 1 - proportion of the read which is aligned, cumulating alignement length on the 2 genomes (here 0.9, meaning 90% of the read has to be aligned)
    # 2 - maximum number of bases inserted between the 2 genomes at recombination point, reflecting non - templated nucleotide additions (here 5)
    # 3 - maximum overlap in the alignement with the 2 genomes at the recombination points, reflects the presence of homology between the 2 genomes at the recombination point (here 20)
    # 4 - minimum alignment length on one genome only. Here the read has to be aligned over at least 16 bp on the genome 1 only and over at least 16 bp on genome 2 only

    mCATblastn_noOverlap$chimeric = chimericReads(mCATblastn_noOverlap, 0.9,  - 5, 20, 16)

    # insert overlap column containing the number of nucleotides shared between the 2 genomes at the recombination point
    mCATblastn_noOverlap$overlap = overlap(mCATblastn_noOverlap)

    # generate data table containing only chimeric reads
    chim_mCATblastn_noOverlap = mCATblastn_noOverlap[chimeric == T]

    # insert column showing the two species involved in the chimeras
    chim_mCATblastn_noOverlap$chimType = ifelse(chim_mCATblastn_noOverlap$Sp>chim_mCATblastn_noOverlap$Sp.s, paste(chim_mCATblastn_noOverlap$Sp, chim_mCATblastn_noOverlap$Sp.s, sep = " - "), paste(chim_mCATblastn_noOverlap$Sp.s, chim_mCATblastn_noOverlap$Sp, sep = " - "))

    # insert column showing which chimeric reads may be PCR duplicates,i.e. which chimeric reads have identical alignement coordinates on both genomes
    chim_mCATblastn_noOverlap[, PCRdup:= paste(qStart - sStart, qEnd - sEnd, subject, qStart.s - sStart.s, qEnd.s - sEnd.s, subject.s)]

    # insert column containing coordinate of the recombination point in one genome
    chim_mCATblastn_noOverlap$chimericPoint = chimericPoint(chim_mCATblastn_noOverlap)

    # insert column containing coordinate of the recombination point in the other genome
    chim_mCATblastn_noOverlap$insertionCoord = insertionCoord(chim_mCATblastn_noOverlap)

    # insert column showing the respective orientation (same or opposite) of the sequences involved in intra - genome chimeras
    chim_mCATblastn_noOverlap[,inv:= sign(sEnd - sStart) != sign(sEnd.s - sStart.s)]

    return(chim_mCATblastn_noOverlap)

}

export_data <- function(df, filename) {
    cat(paste('Exporting data to:', filename, "\n"))
    write.table(df, filename, sep = ',', row.names = FALSE, quote = FALSE)
}

#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


args <- get_args()

blast1 = fread(args$blast_hits_1)
blast2 = fread(args$blast_hits_2)

chim_mCATblastn_noOverlap <- find_chimeras(blast1, blast2)

if ( nrow(chim_mCATblastn_noOverlap) == 0 ) {
    cat("\nNo chimeras found")
} else {
    cat(paste("\nFound ", nrow(chim_mCATblastn_noOverlap), " chimeras"))
    export_data(chim_mCATblastn_noOverlap, args$outfile)
}
