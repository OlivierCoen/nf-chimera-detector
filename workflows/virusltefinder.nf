/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_SRA_IDS                                                     } from '../subworkflows/local/fetch_sra_ids'
include { FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS  as DOWNLOAD_SRA     } from '../subworkflows/local/fastq_download_prefetch_fasterqdump_sratools'
include { POST_PROCESS_SRA                                                  } from '../subworkflows/local/post_process_sra'
include { MMSEQS_WORKFLOW                                                   } from '../subworkflows/local/mmseqs'
include { BLAST_WORKFLOW                                                    } from '../subworkflows/local/blast'
include { MULTIQC_WORKFLOW                                                  } from '../subworkflows/local/multiqc'

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

    def target_db_file = file( params.target_db, checkExists: true )
    ch_target_db = Channel.of([
        [ id: target_db_file.baseName ],
        target_db_file
    ])

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
    POST_PROCESS_SRA.out.single_reads.set { ch_sra_single_reads }

    // ------------------------------------------------------------------------------------
    // PRE-FILTER CANDIDATE QUERY SEQUENCES
    // ------------------------------------------------------------------------------------

    MMSEQS_WORKFLOW (
        ch_sra_single_reads,
        ch_target_db
    )

    PARSE_MMSEQS_OUTPUT (
        MMSEQS_WORKFLOW.out.hits,
        params.max_mmseqs_evalue
    )

    PARSE_MMSEQS_OUTPUT.out.status
        .view { v -> "before fitler" + v}
        .filter { meta, status -> status == "PASS" }
        .view { v -> "after fitler" + v}
        .join( ch_sra_single_reads )
        .map {
            meta, status, sra_read ->
                [ meta, sra_read ]
        }
        .view { v -> "after map" + v}
        .set { ch_filtered_sra_reads }

    // ------------------------------------------------------------------------------------
    // BLAST
    // ------------------------------------------------------------------------------------

    BLAST_WORKFLOW (
        ch_filtered_sra_reads,
        ch_target_db
    )

    BLAST_WORKFLOW.out.hits.view()

    // ------------------------------------------------------------------------------------
    // MULTIQC
    // ------------------------------------------------------------------------------------

    ch_versions
        .mix ( DOWNLOAD_SRA.out.versions )
        .mix ( MMSEQS_WORKFLOW.out.versions )
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
