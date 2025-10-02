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
    ch_reads_fasta
    ch_chimeras_data_mqc


    main:

    ch_data_per_family = Channel.empty()

    // ------------------------------------------------------------------------------------
    // SUMMARY OF CHIMERAS STATISTICS PER FAMILY / SPECIES
    // ------------------------------------------------------------------------------------

    MAKE_SUMMARY (
        ch_chimeras_data_mqc.collect(),
        ['family', 'species']
    )

    // ------------------------------------------------------------------------------------
    // MAKING GENERAL STAT TABLE PER SRR (ID)
    // ------------------------------------------------------------------------------------


    Channel.topic('downloaded_genome_nb_bases')
        .combine( ch_reads_fasta )
        .filter { family, taxid, nb_bases, meta, fasta -> taxid == meta.taxid }
        .map { family, taxid, nb_bases, meta, fasta -> [ meta.id, nb_bases] }
        .set { ch_downloaded_genome_nb_bases }


    Channel.topic('assembled_genome_nb_bases').view{ v -> "assembled_genome_nb_bases " + v}
    Channel.topic('read_fasta_nb_bases').view{ v -> "read_fasta_nb_bases " + v}
    Channel.topic('fastq_sum_len').view{ v -> "fastq_sum_len " + v}

    // ------------------------------------------------------------------------------------
    // COLLECTING NB OF BASES SIZE OF DOWNLOADED FASTQ FILES PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('fastq_sum_len') // family, id, sum_len
        .collectFile(
            name: 'fastq_nb_bases.tsv',
            seed: "family\tdata",
            newLine: true,
            storeDir: "${params.outdir}/sratools/"
        ) {
            item -> "${item[0]}\t${item[2]}"
        }
        .set { ch_fastq_size_file }

    // ------------------------------------------------------------------------------------
    // COLLECTING NB OF BASES SIZE OF PROCESSED READ FASTA FILES PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('read_fasta_nb_bases') // family, id, sum_len
        .collectFile(
            name: 'read_fasta_nb_bases.tsv',
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

    Channel.topic('downloaded_genome_nb_bases') // family, taxid, nb_bases
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

    ch_data_per_family
        .mix ( ch_fastq_size_file )
        .mix ( ch_read_fasta_size_file )
        .mix ( ch_downloaded_genome_size_file )
        .mix ( ch_assembled_genome_size_file )
        .mix ( ch_blast_nb_hits_target_file )
        .mix ( ch_blast_nb_hits_genomes_file )
        .mix ( ch_nb_chimeras_file )
        .set { ch_data_per_family }

    // ------------------------------------------------------------------------------------
    // PIVOTING AND TRANSPOSING DATA FOR MULTIQC
    // ------------------------------------------------------------------------------------

    PREPARE_DATA_PER_FAMILY ( ch_data_per_family )


    emit:
    chimeras_summary                = MAKE_SUMMARY.out.csv
    prepared_data                   = PREPARE_DATA_PER_FAMILY.out.tsv
}

