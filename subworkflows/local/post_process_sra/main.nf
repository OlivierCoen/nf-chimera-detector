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

    // in very rare cases (like SRX18097776), we obtain 3 files from fasterq-dump
    // in suh cases, we keep only the ones having _1/_2 suffix (especially for convenience)
    ch_branched_sra_reads.paired
        .map {
            meta, reads ->
                if ( reads.size() == 3 ) {
                    kept_reads = reads.findAll { it ==~ /.*_[12]\.f(ast)?q\.gz$/  }
                    [ meta, kept_reads ]
                } else {
                    [ meta, reads ]
                }
        }
        .set { ch_sra_paired_reads }

    def interleave = false
    BBMAP_BBMERGE (
        ch_sra_paired_reads,
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

