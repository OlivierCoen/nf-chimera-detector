include { BLAST_MAKEBLASTDB                     } from '../../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN                          } from '../../../modules/nf-core/blast/blastn'


workflow BLAST_WORKFLOW {

    take:
    ch_sra_reads
    ch_target_db

    main:

    ch_versions = Channel.empty()

    BLAST_MAKEBLASTDB ( ch_target_db )
    ch_sra_reads.view()
    BLAST_BLASTN (
        ch_sra_reads,
        BLAST_MAKEBLASTDB.out.db
    )

    ch_versions
        .mix ( BLAST_MAKEBLASTDB.out.versions )
        .mix ( BLAST_BLASTN.out.versions )
        .set { ch_versions }

    emit:
    hits                            = BLAST_BLASTN.out.txt
    versions                        = ch_versions

}

