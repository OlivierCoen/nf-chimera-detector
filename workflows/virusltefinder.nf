/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_SRA_IDS                                           } from '../subworkflows/local/fetch_sra_ids'
include { MULTIQC_WORKFLOW                                        } from '../subworkflows/local/multiqc'
include { MMSEQS_WORKFLOW                                         } from '../subworkflows/local/mmseqs'
include { FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS            } from '../subworkflows/local/fastq_download_prefetch_fasterqdump_sratools'



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

    // ------------------------------------------------------------------------------------
    // GETTING LIST OF SRA IDS
    // ------------------------------------------------------------------------------------

    FETCH_SRA_IDS ( ch_families )

    // ------------------------------------------------------------------------------------
    // DOWNLOAD ALL SRA DATA
    // ------------------------------------------------------------------------------------

    FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS (
        FETCH_SRA_IDS.out.sra_ids,
        params.dbgap_key ? file(params.dbgap_key, checkIfExists: true) : []
    )

    // ------------------------------------------------------------------------------------
    // PRE-FILTER CANDIDATE QUERY SEQUENCES
    // ------------------------------------------------------------------------------------

    ch_target_db = Channel.of([
        [ id: "db" ],
        file( params.target_db, checkExists: true )
    ])

    MMSEQS_WORKFLOW (
        FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS.out.reads,
        ch_target_db
    )

    // ------------------------------------------------------------------------------------
    // MULTIQC
    // ------------------------------------------------------------------------------------

    ch_versions
        .mix ( FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS.out.versions )
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
