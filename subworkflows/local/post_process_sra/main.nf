include { BBMAP_BBMERGE                         } from '../../../modules/nf-core/bbmap/bbmerge'
include { SEQKIT_FQ2FA                          } from '../../../modules/local/seqkit/fq2fa'


workflow POST_PROCESS_SRA {

    take:
    ch_sra_reads

    main:

    ch_versions = Channel.empty()

    ch_sra_reads
        .branch {
            meta, reads ->
                single: reads instanceof Path
                paired: reads instanceof List
        }
        .set { ch_branched_sra_reads }

    def interleave = false
    BBMAP_BBMERGE (
        ch_branched_sra_reads.paired,
        interleave
    )

    ch_branched_sra_reads.single
        .mix ( BBMAP_BBMERGE.out.merged )
        .set { ch_sra_single_reads }

    SEQKIT_FQ2FA ( ch_sra_single_reads )

    ch_versions
        .mix ( BBMAP_BBMERGE.out.versions )
        .set { ch_versions }

    emit:
    single_reads                    = SEQKIT_FQ2FA.out.fasta
    versions                        = ch_versions

}

