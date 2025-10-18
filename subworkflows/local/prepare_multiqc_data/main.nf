include { PREPARE_DATA_PER_FAMILY                           } from '../../../modules/local/prepare_data_per_family'
include { MAKE_SUMMARY                                      } from '../../../modules/local/make_summary'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def getBlastDb( process_name ) {
    return process_name.tokenize(':')[1].tokenize('_')[2]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow PREPARE_MULTIQC_DATA {

    take:
    ch_chimeras_table
    ch_fastq_stats
    ch_reads_fasta
    ch_species_taxids


    main:

    ch_data_per_family = Channel.empty()

    // ------------------------------------------------------------------------------------
    // PREPARING FASTQ STATISTICS FOR THE GENERALSTATS TABLE
    // ------------------------------------------------------------------------------------

    ch_fastq_stats
        .collectFile(
            name: 'fastq_stats.tsv',
            seed: "srr_id\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\tsum_gap\tN50\tN50_num\tQ20_pct\tQ30_pct\tAvgQual\tGC_pct\tsum_n",
            newLine: true,
            storeDir: "${params.outdir}/multiqc/"
        ) {
            item -> "${item.id}\t${item.num_seqs}\t${item.sum_len}\t${item.min_len}\t${item.avg_len}\t${item.max_len}\t${item.Q1}\t${item.Q2}\t${item.Q3}\t${item.sum_gap}\t${item.N50}\t${item.N50_num}\t${item['Q20(%)']}\t${item['Q30(%)']}\t${item.AvgQual}\t${item['GC(%)']}\t${item.sum_n}"
        }
        .set { ch_fastq_stats_file }

    // ------------------------------------------------------------------------------------
    // GETTING NB OF SPECIES PER FAMILY
    // ------------------------------------------------------------------------------------

    ch_species_taxids
        .map { meta, taxid -> [ meta.family, taxid ] }
        .groupTuple()
        .map { meta, taxids -> [ meta, taxids.size() ] }
        .collectFile(
            name: 'nb_species_per_family.tsv',
            seed: "family\tnb_species",
            newLine: true,
            storeDir: "${params.outdir}/species_taxids/"
        ) {
            item -> "${item[0]}\t${item[1]}"
        }
        .set { ch_nb_species_per_family_file }

    // ------------------------------------------------------------------------------------
    // GETTING NB SRRS PER FAMILY
    // ------------------------------------------------------------------------------------

    ch_reads_fasta
        .map { meta, file -> [ meta.family, file ] }
        .groupTuple()
        .map { meta, files -> [ meta, files.size() ] }
        .collectFile(
            name: 'nb_srrs_per_family.tsv',
            seed: "family\tnb_srrs",
            newLine: true,
            storeDir: "${params.outdir}/sratools/"
        ) {
            item -> "${item[0]}\t${item[1]}"
        }
        .set { ch_nb_srrs_per_family_file }


    // ------------------------------------------------------------------------------------
    // SUMMARY OF CHIMERAS STATISTICS PER FAMILY / SPECIES
    // ------------------------------------------------------------------------------------

    MAKE_SUMMARY (
        ch_chimeras_table,
        ['family', 'species']
    )

    // ------------------------------------------------------------------------------------
    // MAKING GENERAL STAT TABLE PER SRR (ID)
    // ------------------------------------------------------------------------------------

    // associating metadata of downloaded genomes to all possible SRR IDs
    Channel.topic('dl_genome_metadata')
        .combine( ch_reads_fasta )
        .filter { family, taxid, genome_file, nb_bases, meta, fasta -> taxid == meta.taxid }
        .map { family, taxid, genome_file, nb_bases, meta, fasta ->  [ meta.id, genome_file.baseName, nb_bases ] }
        .set { ch_dl_genome_len }

    // performing multiple left joins (remainder is true) in a row, with id (SRR ID) as matching key
    ch_reads_fasta
        .map { meta, fasta -> [ meta.id, meta ] } // putting SRR id separately in order to join easily with other channels
        .join( // joining with assembly metadata
            Channel.topic('asm_genome_len').map { family, id, nb -> [ id, [asm_genome_len: nb] ] },
            remainder: true
        )
        .join( // joining with read fasta metadata
            Channel.topic('read_fasta_len').map { family, id, nb -> [ id, [read_fasta_len: nb] ] },
            remainder: true
        )
        .join ( // joining with download genome metadata
            ch_dl_genome_len.map { id, genome_name, nb -> [ id, [dl_genome_name: genome_name, dl_genome_len: nb] ] },
            remainder: true
        )
        .join( // joining with chimeras metadata
            Channel.topic('nb_chimeras').map { family, id, nb -> [ id, [nb_chimeras: nb] ] },
            remainder: true
        )
        .map { // removing SRR ID and cleaning data
            meta ->
                // removes first element: taxid (1 does not correspond to the position but to the nb of elements removes from the beginning)
                def cleaned = meta.drop(1)
                // drop nulls (only for genomes: downloaded XOR assembled)
                def not_nulls = cleaned.findAll { it != null }
                // merge list of maps into one map
                def merged  = not_nulls.collectEntries { it }
                merged
        }
        .map { // computing coverage
            meta ->
                def genome_length = meta.dl_genome_len ?: meta.asm_genome_len
                def coverage = meta.read_fasta_len.toFloat() / genome_length.toFloat()
                meta + [coverage: coverage]
        }
        .collectFile(
            name: 'srr_metadata.tsv',
            seed: "srr_id\tfamily\ttaxid\ttaxon_name\tsra_id\tcoverage\tread_fasta_len\tdl_genome_name\tdl_genome_len\tasm_genome_len\tnb_chimeras", // header of TSV file
            newLine: true,
            storeDir: "${params.outdir}/multiqc/"
        ) {
            item -> "${item.id}\t${item.family}\ttxid${item.taxid}\t${item.taxon_name}\t${item.sra_id}\t${item.coverage}\t${item.read_fasta_len}\t${item.dl_genome_name}\t${item.dl_genome_len}\t${item.asm_genome_len}\t${item.nb_chimeras}"
        }
        .set { ch_srr_metadata_file }

    // ------------------------------------------------------------------------------------
    // COLLECTING NB OF BASES SIZE OF PROCESSED READ FASTA FILES PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('read_fasta_len') // family, id, sum_len
        .collectFile(
            name: 'read_fasta_len.tsv',
            seed: "family\tdata",
            newLine: true,
            storeDir: "${params.outdir}/sratools/"
        ) {
            item -> "${item[0]}\t${item[2]}"
        }
        .set { ch_read_fasta_size_file }

    // ------------------------------------------------------------------------------------
    // COLLECTING NB OF BASES OF DOWNLOADED GENOMES PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('dl_genome_len') // family, taxid, nb_bases
        .collectFile(
            name: 'downloaded_genome_size.tsv',
            seed: "family\tdata",
            newLine: true,
            storeDir: "${params.outdir}/assemblies/"
        ) {
            item -> "${item[0]}\t${item[2]}"
        }
        .set { ch_downloaded_genome_size_file }

    // ------------------------------------------------------------------------------------
    // COLLECTING NB OF BASES OF ASSEMBLED GENOMES PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('assembled_genome_size') // family, id, nb_bases
        .collectFile(
            name: 'assembled_genome_size.tsv',
            seed: "family\tdata",
            newLine: true,
            storeDir: "${params.outdir}/assemblies/"
        ) {
            item -> "${item[0]}\t${item[2]}"
        }
        .set { ch_assembled_genome_size_file }

    // ------------------------------------------------------------------------------------
    // COLLECTING NB OF BLAST HITS AGAINST TARGET / GENOME PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('blast_nb_hits')
        .branch {
            process, family, id, nb_hits ->
                target: getBlastDb(process) == 'TARGET'
                    [ family, id, nb_hits ]
                genomes: getBlastDb(process) == 'GENOMES'
                   [ family, id, nb_hits ]
        }
        .set { ch_blast_nb_hits }

    ch_blast_nb_hits.target // family, id, nb_hits
        .collectFile(
            name: 'blast_nb_hits_target.tsv',
            seed: "family\tdata",
            newLine: true,
            storeDir: "${params.outdir}/blastn/"
        ) {
            item -> "${item[0]}\t${item[2]}"
        }
        .set { ch_blast_nb_hits_target_file }


    ch_blast_nb_hits.genomes // family, id, nb_hits
        .collectFile(
            name: 'blast_nb_hits_genomes.tsv',
            seed: "family\tdata",
            newLine: true,
            storeDir: "${params.outdir}/blastn/"
        ) {
            item -> "${item[0]}\t${item[2]}"
        }
        .set { ch_blast_nb_hits_genomes_file }

    // ------------------------------------------------------------------------------------
    // COLLECTING NB OF CHIMERAS PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('nb_chimeras')  // family, id, nb_chimeras
        .collectFile(
            name: 'nb_chimeras.tsv',
            seed: "family\tdata",
            newLine: true,
            storeDir: "${params.outdir}/chimeras/"
        ) {
            item -> "${item[0]}\t${item[2]}"
        }
        .set { ch_nb_chimeras_file }

    // ------------------------------------------------------------------------------------
    // PIVOTING AND TRANSPOSING DATA FOR MULTIQC
    // ------------------------------------------------------------------------------------

    ch_data_per_family
        .mix ( ch_read_fasta_size_file )
        .mix ( ch_downloaded_genome_size_file )
        .mix ( ch_assembled_genome_size_file )
        .mix ( ch_blast_nb_hits_target_file )
        .mix ( ch_blast_nb_hits_genomes_file )
        .mix ( ch_nb_chimeras_file )
        .set { ch_data_per_family }

    PREPARE_DATA_PER_FAMILY ( ch_data_per_family )


    emit:
    srr_metadata                    = ch_srr_metadata_file
    fastq_stats                     = ch_fastq_stats_file
    chimeras_summary                = MAKE_SUMMARY.out.csv
    nb_species_per_family           = ch_nb_species_per_family_file
    nb_srrs_per_family              = ch_nb_srrs_per_family_file
    data_per_family                 = PREPARE_DATA_PER_FAMILY.out.tsv
}
