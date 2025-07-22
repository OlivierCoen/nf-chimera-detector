include { BLAST_MAKEBLASTDB                     } from '../../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN                          } from '../../../modules/local/blast/blastn'


workflow BLAST_WORKFLOW {

    take:
    ch_sra_reads
    ch_target_db

    main:

    ch_versions = Channel.empty()

    BLAST_MAKEBLASTDB ( ch_target_db )

     // making all combinations of reads + target db
    ch_sra_reads
        .combine( BLAST_MAKEBLASTDB.out.db )
        .map { meta, reads, meta2, db ->  [ meta, reads, db ] }
        .set { blastn_input }

    BLAST_BLASTN ( blastn_input )

    ch_versions
        .mix ( BLAST_MAKEBLASTDB.out.versions )
        .set { ch_versions }

    emit:
    hits                            = BLAST_BLASTN.out.txt
    versions                        = ch_versions

}

