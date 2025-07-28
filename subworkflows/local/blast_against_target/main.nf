include { BLAST_MAKEBLASTDB as MAKEBLASTDB             } from '../../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN as BLASTN                       } from '../../../modules/local/blast/blastn'


workflow BLAST_AGAINST_TARGET {

    take:
    ch_sra_reads
    ch_target_db

    main:

    ch_versions = Channel.empty()

    MAKEBLASTDB ( ch_target_db )

     // making all combinations of reads + target db
    ch_sra_reads
        .combine( MAKEBLASTDB.out.db )
        .map { meta, reads, meta2, db ->  [ meta, reads, db ] }
        .set { blastn_input }

    BLASTN ( blastn_input )

    ch_versions
        .mix ( MAKEBLASTDB.out.versions )
        .set { ch_versions }

    emit:
    hits                            = BLASTN.out.txt
    versions                        = ch_versions

}

