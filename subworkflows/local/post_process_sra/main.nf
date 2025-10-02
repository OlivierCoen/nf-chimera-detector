include { SEQKIT_PAIR                           } from '../../../modules/nf-core/seqkit/pair'
include { BBMAP_BBMERGE as BBMERGE              } from '../../../modules/nf-core/bbmap/bbmerge'
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
    // in such cases, we keep only the ones having _1/_2 suffix (especially for convenience)
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
    // FILTER OUT UNPAIRED READS (NECESSARY FOR BBMERGE)
    // ------------------------------------------------------------------------------------

    SEQKIT_PAIR ( ch_sra_paired_reads )

    // ------------------------------------------------------------------------------------
    // MERGE OVERLAPPING PAIRED READS INTO A SINGLE FASTQ FILE
    // ------------------------------------------------------------------------------------

    def interleave = false
    BBMERGE (
        SEQKIT_PAIR.out.reads,
        interleave
    )

    // putting together merged and unmerged reads
    BBMERGE.out.merged
        .join( BBMERGE.out.unmerged )
        .map{
            meta, merged, unmerged ->
                def reads = [ merged, unmerged ]
                [ meta, reads.flatten() ]
        }
        .set { ch_processed_paired_reads }

    // putting together single reads and processed paired reads
    ch_branched_sra_reads.single
        .mix ( ch_processed_paired_reads )
        .set { ch_reads }

    // ------------------------------------------------------------------------------------
    // FASTQ TO FASTA
    // ------------------------------------------------------------------------------------

    SEQKIT_FQ2FA ( ch_reads )


    ch_versions
        .mix ( SEQKIT_PAIR.out.versions )
        .mix ( BBMERGE.out.versions )
        .set { ch_versions }

    emit:
    merged_reads_fasta              = SEQKIT_FQ2FA.out.fasta
    versions                        = ch_versions

}

