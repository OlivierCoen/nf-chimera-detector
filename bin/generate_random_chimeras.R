###########################################################################################################################
########### Generate random artificial chimeric reads between regions of the same genome or between two genomes ###########
###########################################################################################################################

#this is to generate a null expectation for the lengths of microhomologies in chimeric reads

require(Biostrings)
require(stringi)

# import genome 1
Genome1 = readDNAStringSet("Genome1.fasta")

Genome1 = stri_flatten(Genome1)   #Caution. Total contig length should not exceed 2^31 bases.

# import genome 2
Genome2 = readDNAStringSet("Genome2.fasta")

Genome2 = stri_flatten(Genome2)    #Total contig length should not exceed 2^31 bases.

# Generate random start coordinate of 1st region of 100000 chimeric reads
startRegion1 = sample(1:(nchar(Genome1) - 150), 100000, T)

# Generate random start coordinate of 2nd region of 100000 chimeric reads
startRegion2 = sample(1:(nchar(Genome2) - 150), 100000, T)   # To generate random chimeras within the same genome, use Genome1 instead of Genome2 here

# Generate random length of 1st region of the chimeric reads - here set to be between 28 bp (minimum length of alignment returned by blastn) and 150 bp (length of the reads)
length1 = sample(28:(150 - 28), 100000, T)

# Generate random length of 2nd region of the chimeric reads
length2 = 150 - length1

# Extract 100000 1st region of the chimeric reads from genome 1
part1 = stri_sub(Genome1, from = startRegion1, length = length1)

# Extract 100000 2nd region of the chimeric reads from genome 2 (or genome 1 when intra - genome chimeras are generated)
part2 = stri_sub(Genome2, from = startRegion2, length = length2)  # To generate random chimeras within the same genome, use Genome1 instead of Genome2 here

# Generate 100000/2 random numbers "N1" of 1st regions to reverse
toReverse1 = sample(1:100000, 100000/2)

# Generate 100000/2 random numbers of 2nd regions to reverse
toReverse2 = sample(1:100000, 100000/2)

# Reverse 100000/2 1st regions
part1[toReverse1] = revCom(part1[toReverse1])

# Reverse 100000/2 2nd regions
part2[toReverse2] = revCom(part2[toReverse2])

# Paste 100000 1st and 100000 2nd regions together
RandomChimericReads = stri_c(part1, part2)

# Name random chimeric reads "read_1" to "read_100000"
names(RandomChimericReads) = paste("read",1:length(RandomChimericReads), sep = "_")

# Write fasta file containing random chimeric reads in working directory
writeXStringSet(DNAStringSet(RandomChimericReads), "RandomChimericReads.fasta")

# The fasta file can then be used as query to blast the various genomes between which recombination is studied. The output of this blast can then be analysed using the above
# script to characterize the features of chimeric reads that are expected to occurr if chimeras are produced randomly.
