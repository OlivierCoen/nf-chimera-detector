include { PREPARE_DATA_PER_FAMILY as PREPARE_FASTQ_SIZE                     } from '../../../modules/local/prepare_data_per_family'
include { PREPARE_DATA_PER_FAMILY as PREPARE_DOWNLOADED_GENOME_SIZE         } from '../../../modules/local/prepare_data_per_family'
include { PREPARE_DATA_PER_FAMILY as PREPARE_ASSEMBLED_GENOME_SIZE          } from '../../../modules/local/prepare_data_per_family'
include { PREPARE_DATA_PER_FAMILY as PREPARE_BLAST_HIT_TARGET               } from '../../../modules/local/prepare_data_per_family'
include { PREPARE_DATA_PER_FAMILY as PREPARE_BLAST_HIT_GENOMES              } from '../../../modules/local/prepare_data_per_family'
include { PREPARE_DATA_PER_FAMILY as PREPARE_NB_CHIMERAS                    } from '../../../modules/local/prepare_data_per_family'

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


    main:

    // ------------------------------------------------------------------------------------
    // DISTRIBUTION OF SIZE OF DOWNLOADED FASTQ FILES (IN NB OF BASES) PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('fastq_size')
        .collectFile(
            name: 'fastq_size.csv',
            seed: "family,data",
            newLine: true,
            storeDir: "${params.outdir}/sratools/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_fastq_size_file }

    PREPARE_FASTQ_SIZE ( ch_fastq_size_file )

    // ------------------------------------------------------------------------------------
    // DISTRIBUTION OF SIZE OF DOWNLOADED GENOMES PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('downloaded_genome_size')
        .collectFile(
            name: 'downloaded_genome_size.csv',
            seed: "family,data",
            newLine: true,
            storeDir: "${params.outdir}/assemblies/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_downloaded_genome_size_file }

    PREPARE_DOWNLOADED_GENOME_SIZE ( ch_downloaded_genome_size_file )

    // ------------------------------------------------------------------------------------
    // DISTRIBUTION OF SIZE OF ASSEMBLED GENOMES PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('assembled_genome_size')
        .collectFile(
            name: 'assembled_genome_size.csv',
            seed: "family,data",
            newLine: true,
            storeDir: "${params.outdir}/assemblies/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_assembled_genome_size_file }

    PREPARE_ASSEMBLED_GENOME_SIZE ( ch_assembled_genome_size_file )

    // ------------------------------------------------------------------------------------
    // DISTRIBUTION OF NB OF BLAST HITS AGAINST TARGET PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('blast_nb_hits')
        .branch {
            process, family, nb_hits ->
                target: getBlastDb(process) == 'TARGET'
                    [ family, nb_hits ]
                genomes: getBlastDb(process) == 'GENOMES'
                   [ family, nb_hits ]
        }
        .set { ch_blast_nb_hits }

    ch_blast_nb_hits.target
        .collectFile(
            name: 'blast_nb_hits_target.csv',
            seed: "family,data",
            newLine: true,
            storeDir: "${params.outdir}/blastn/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_blast_nb_hits_target_file }

    PREPARE_BLAST_HIT_TARGET ( ch_blast_nb_hits_target_file )

    // ------------------------------------------------------------------------------------
    // DISTRIBUTION OF NB OF BLAST HITS AGAINST GENOMES PER FAMILY
    // ------------------------------------------------------------------------------------

    ch_blast_nb_hits.genomes
        .collectFile(
            name: 'blast_nb_hits_genomes.csv',
            seed: "family,data",
            newLine: true,
            storeDir: "${params.outdir}/blastn/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_blast_nb_hits_genomes_file }

    PREPARE_BLAST_HIT_GENOMES ( ch_blast_nb_hits_genomes_file )

    // ------------------------------------------------------------------------------------
    // DISTRIBUTION OF NB OF CHIMERAS PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('nb_chimeras')
        .collectFile(
            name: 'nb_chimeras.csv',
            seed: "family,data",
            newLine: true,
            storeDir: "${params.outdir}/chimeras/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_nb_chimeras_file }

    PREPARE_NB_CHIMERAS ( ch_nb_chimeras_file )


    emit:
    fastq_sizes                     = PREPARE_FASTQ_SIZE.out.csv
    downloaded_genome_sizes         = PREPARE_DOWNLOADED_GENOME_SIZE.out.csv
    assembled_genome_sizes          = PREPARE_ASSEMBLED_GENOME_SIZE.out.csv
    nb_blast_hits_target            = PREPARE_BLAST_HIT_TARGET.out.csv
    nb_blast_hits_genomes           = PREPARE_BLAST_HIT_GENOMES.out.csv
    nb_chimeras                     = PREPARE_NB_CHIMERAS.out.csv
}

