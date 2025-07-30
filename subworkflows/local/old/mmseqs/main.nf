include { MMSEQS_CREATE_DB_WITH_INDEX                } from '../../../modules/local/mmseqs/create_db_with_index'
include { MMSEQS_EASYSEARCH                          } from '../../../modules/local/mmseqs/easysearch'


workflow MMSEQS_WORKFLOW {

    take:
    ch_sra_reads
    ch_target_db

    main:

    ch_versions = Channel.empty()

    MMSEQS_CREATE_DB_WITH_INDEX ( ch_target_db )

    // making all combinations of reads + target db
    ch_sra_reads
        .combine( MMSEQS_CREATE_DB_WITH_INDEX.out.db )
        .map { meta, reads, meta2, db ->  [ meta, reads, db ] }
        .set { mmseqs_easysearch_input }

    MMSEQS_EASYSEARCH ( mmseqs_easysearch_input )

    emit:
    hits                            = MMSEQS_EASYSEARCH.out.tsv
    versions                        = ch_versions

}

