/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_SRA_IDS                                                             } from '../subworkflows/local/fetch_sra_ids'
include { DOWNLOAD_SRA                                                              } from '../subworkflows/local/download_sra'
include { ASSEMBLY                                                                  } from '../subworkflows/local/assembly'
include { POST_PROCESS_SRA                                                          } from '../subworkflows/local/post_process_sra'
include { MMSEQS_WORKFLOW                                                           } from '../subworkflows/local/mmseqs'
include { BLAST_AGAINST_TARGET                                                      } from '../subworkflows/local/blast_against_target'
include { BLAST_AGAINST_ASSEMBLIES                                                  } from '../subworkflows/local/blast_against_assemblies'
include { MULTIQC_WORKFLOW                                                          } from '../subworkflows/local/multiqc'

include { NCBI_ASSEMBLY_STATS                                                       } from '../modules/local/ncbi_assembly_stats'
include { PARSE_MMSEQS_OUTPUT                                                       } from '../modules/local/parse_mmseqs_output'

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

    def target_fasta_file = file( params.target_db, checkExists: true )
    ch_target_db = Channel.of([
        [ id: target_fasta_file.baseName ],
        target_fasta_file
    ])

    NCBI_ASSEMBLY_STATS ( ch_families )

    NCBI_ASSEMBLY_STATS.out.mean_lengths
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
    DOWNLOAD_SRA.out.reads.set { ch_sra_reads }

    // ------------------------------------------------------------------------------------
    // DOWNLOAD GENOME ASSEMBLY FROM NCBI IF AVAILABLE OR MAKE FROM SCRATCH OTHERWISE
    // ------------------------------------------------------------------------------------

    ASSEMBLY ( ch_sra_reads )
    ASSEMBLY.out.assemblies.set { ch_assemblies }

    // ------------------------------------------------------------------------------------
    // COMBINE PAIRED READS (IF NECESSARY) AND CONVERT FASTQ TO FASTA
    // ------------------------------------------------------------------------------------

    POST_PROCESS_SRA ( ch_sra_reads )
    POST_PROCESS_SRA.out.single_reads.set { ch_sra_reads }

    // ------------------------------------------------------------------------------------
    // FILTER OUT SRA READS NOT SHOWING ANY HIT AGAINST TARGET DB
    // ------------------------------------------------------------------------------------

    MMSEQS_WORKFLOW (
        ch_sra_reads,
        ch_target_db
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

    PARSE_MMSEQS_OUTPUT.out.status
        .filter { meta, status -> status == "FAIL" }
        .view { meta, _ -> "Failed MMSeqs filtering: Family: ${meta.family} | Taxid: ${meta.taxid} | SRA ID: ${meta.id}" }

    // ------------------------------------------------------------------------------------
    // BLAST AGAINST TARGET
    // ------------------------------------------------------------------------------------

    BLAST_AGAINST_TARGET (
        ch_sra_reads,
        ch_target_db
    )
    BLAST_AGAINST_TARGET.out.hits.set { ch_blast_target_hits }

    // ------------------------------------------------------------------------------------
    // BLAST AGAINST ASSEMBLIES
    // ------------------------------------------------------------------------------------

    BLAST_AGAINST_ASSEMBLIES (
        ch_blast_target_hits,
        ch_assemblies
    )
    BLAST_AGAINST_ASSEMBLIES.out.hits.set { ch_blast_assembly_hits }

    // ------------------------------------------------------------------------------------
    // MULTIQC
    // ------------------------------------------------------------------------------------

    ch_versions
        .mix ( DOWNLOAD_SRA.out.versions )
        .mix ( BLAST_AGAINST_ASSEMBLIES.out.versions )
        .mix ( MMSEQS_WORKFLOW.out.versions )
        .set { ch_versions }

    //MULTIQC_WORKFLOW ( ch_versions )
    ch_a = Channel.empty()
    emit:
    multiqc_report = ch_a
    //multiqc_report = MULTIQC_WORKFLOW.out.multiqc_report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
