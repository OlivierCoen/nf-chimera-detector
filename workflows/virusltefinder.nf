/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_SRA_IDS                                                     } from '../subworkflows/local/fetch_sra_ids'
include { FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS  as DOWNLOAD_SRA     } from '../subworkflows/local/fastq_download_prefetch_fasterqdump_sratools'
include { POST_PROCESS_SRA                                                  } from '../subworkflows/local/post_process_sra'
include { MMSEQS_WORKFLOW                                                   } from '../subworkflows/local/mmseqs'
include { BLAST_WORKFLOW as BLAST_AGAINST_TE                                } from '../subworkflows/local/blast'
include { MULTIQC_WORKFLOW                                                  } from '../subworkflows/local/multiqc'

include { GET_MEAN_ASSEMBLY_LENGTH                                          } from '../modules/local/get_mean_assembly_length'
include { PARSE_MMSEQS_OUTPUT                                               } from '../modules/local/parse_mmseqs_output'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VIRUSLTEFINDER {

    take:
    ch_families

    main:

    ch_versions = Channel.empty()

    def te_fasta_file = file( params.te_fasta, checkExists: true )
    ch_te_db = Channel.of([
        [ id: te_fasta_file.baseName ],
        te_fasta_file
    ])

    GET_MEAN_ASSEMBLY_LENGTH ( ch_families )

    GET_MEAN_ASSEMBLY_LENGTH.out.families
        .map {
            family, mean_assembly_length ->
                def meta = [ family: family, mean_assembly_length: mean_assembly_length ]
                [ meta, family ]
        }
        .set { ch_families }


    // ------------------------------------------------------------------------------------
    // GETTING LIST OF SRA IDS
    // ------------------------------------------------------------------------------------

    FETCH_SRA_IDS ( ch_families )

    // ------------------------------------------------------------------------------------
    // DOWNLOAD ALL SRA DATA
    // ------------------------------------------------------------------------------------

    DOWNLOAD_SRA ( FETCH_SRA_IDS.out.sra_ids )

    // ------------------------------------------------------------------------------------
    // COMBINE PAIRED READS AND CONVERT TO FASTA
    // ------------------------------------------------------------------------------------

    POST_PROCESS_SRA ( DOWNLOAD_SRA.out.reads )
    POST_PROCESS_SRA.out.single_reads.set { ch_sra_reads }

    // ------------------------------------------------------------------------------------
    // PRE-FILTER CANDIDATE QUERY SEQUENCES
    // ------------------------------------------------------------------------------------

    if ( !params.skip_mmseqs_prefiltering ) {

        MMSEQS_WORKFLOW (
            ch_sra_reads,
            ch_te_db
        )

        PARSE_MMSEQS_OUTPUT (
            MMSEQS_WORKFLOW.out.hits,
            params.max_mmseqs_evalue
        )

        PARSE_MMSEQS_OUTPUT.out.status
            .filter { meta, status -> status == "PASS" }
            .join( ch_sra_reads )
            .map {
                meta, status, sra_read ->
                    [ meta, sra_read ]
            }
            .set { ch_sra_reads }

        ch_versions = ch_versions.mix ( MMSEQS_WORKFLOW.out.versions )

    }

    // ------------------------------------------------------------------------------------
    // BLAST AGAINST TE representative sequences
    // ------------------------------------------------------------------------------------

    BLAST_AGAINST_TE (
        ch_sra_reads,
        ch_te_db
    )
    BLAST_AGAINST_TE.out.hits.set { ch_blast_te_hits }

    /*
    // filtering reads that had a hit
    ch_sra_reads
        .join ( ch_blast_te_hits )
        .map { meta, reads, txt -> [ meta, reads ]}


    ch_blast_te_hits.view()

    def target_genome_fasta_file = file( params.genomes_fasta, checkExists: true )
    ch_target_genome_db = Channel.of([
        [ id: target_genome_fasta_file.baseName ],
        target_genome_fasta_file
    ])
    */

    // ------------------------------------------------------------------------------------
    // MULTIQC
    // ------------------------------------------------------------------------------------

    ch_versions
        .mix ( DOWNLOAD_SRA.out.versions )
        .mix ( BLAST_AGAINST_TE.out.versions )
        .set { ch_versions }

    MULTIQC_WORKFLOW ( ch_versions )

    emit:
    multiqc_report = MULTIQC_WORKFLOW.out.multiqc_report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
