include { SEQKIT_PAIR                           } from '../../../modules/nf-core/seqkit/pair'
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

    // ------------------------------------------------------------------------------------
    // FILTER OUT UNPAIRED READS
    // ------------------------------------------------------------------------------------

    SEQKIT_PAIR ( ch_sra_paired_reads )

    // ------------------------------------------------------------------------------------
    // MERGE PAIRED READS INTO A SINGLE FASTQ FILE
    // ------------------------------------------------------------------------------------

    def interleave = false
    BBMAP_BBMERGE (
        SEQKIT_PAIR.out.reads,
        interleave
    )

    ch_branched_sra_reads.single
        .mix ( BBMAP_BBMERGE.out.merged )
        .set { ch_sra_single_reads }

    // ------------------------------------------------------------------------------------
    // FASTQ TO FASTA
    // ------------------------------------------------------------------------------------

    SEQKIT_FQ2FA ( ch_sra_single_reads )

    ch_versions
        .mix ( SEQKIT_PAIR.out.versions )
        .mix ( BBMAP_BBMERGE.out.versions )
        .set { ch_versions }

    emit:
    merged_reads_fasta              = SEQKIT_FQ2FA.out.fasta
    versions                        = ch_versions

}

