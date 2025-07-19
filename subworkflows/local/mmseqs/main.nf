include { MMSEQS_CREATEDB                            } from '../../../modules/nf-core/mmseqs/createdb'
include { MMSEQS_CREATEINDEX                         } from '../../../modules/nf-core/mmseqs/createindex'
include { MMSEQS_EASYSEARCH                          } from '../../../modules/nf-core/mmseqs/easysearch'


workflow MMSEQS_WORKFLOW {

    take:
    ch_sra_reads
    ch_target_db

    main:

    ch_versions = Channel.empty()

    MMSEQS_CREATEDB ( ch_target_db )

    MMSEQS_CREATEINDEX ( MMSEQS_CREATEDB.out.db )

    MMSEQS_EASYSEARCH (
        ch_sra_reads,
        ch_target_db
    )

    ch_versions
        .mix ( MMSEQS_CREATEDB.out.versions )
        .mix ( MMSEQS_CREATEINDEX.out.versions )
        .mix ( MMSEQS_EASYSEARCH.out.versions )
        .set { ch_versions }

    emit:
    hits                            = MMSEQS_EASYSEARCH.out.tsv
    versions                        = ch_versions

}

