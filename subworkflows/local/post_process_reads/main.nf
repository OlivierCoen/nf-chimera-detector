include { SEQKIT_PAIR                           } from '../../../modules/local/seqkit/pair'
include { BBMAP_BBMERGE as BBMERGE              } from '../../../modules/nf-core/bbmap/bbmerge'
include { FLASH                                 } from '../../../modules/local/flash'
include { SEQKIT_FQ2FA                          } from '../../../modules/local/seqkit/fq2fa'


workflow POST_PROCESS_READS {

    take:
    ch_reads

    main:

    ch_versions = channel.empty()

    ch_branched_reads = ch_reads
                            .branch {
                                meta, reads ->
                                    single: reads instanceof Path
                                    paired: reads instanceof List
                            }

    // in very rare cases (like SRX18097776), we obtain 3 files from fasterq-dump
    // in such cases, we keep only the ones having _1/_2 suffix (especially for convenience)
    ch_paired_reads = ch_branched_reads.paired
                        .map {
                            meta, reads ->
                                if ( reads.size() == 3 ) {
                                    kept_reads = reads.findAll { it ==~ /.*_[12]\.f(ast)?q\.gz$/  }
                                    //log.warn "SRR " + meta.id + " has 3 downloaded files!"
                                    [ meta, kept_reads ]
                                } else {
                                    [ meta, reads ]
                                }
                        }


    // ------------------------------------------------------------------------------------
    // FILTER OUT UNPAIRED READS (NECESSARY FOR BBMERGE)
    // ------------------------------------------------------------------------------------

    SEQKIT_PAIR ( ch_paired_reads )
    ch_properly_paired_reads = SEQKIT_PAIR.out.reads

    // ------------------------------------------------------------------------------------
    // MERGE OVERLAPPING PAIRED READS INTO A SINGLE FASTQ FILE
    // ------------------------------------------------------------------------------------

    if ( params.read_merger == "flash" ) {

        FLASH ( ch_properly_paired_reads )

        // putting together merged and unmerged reads
        ch_processed_paired_reads = FLASH.out.merged
                                        .join( FLASH.out.notcombined )
                                        .map{
                                            meta, merged, unmerged ->
                                                def reads = [ merged, unmerged ]
                                                [ meta, reads.flatten() ]
                                        }

    } else if ( params.read_merger == "bbmerge" ) {

        def interleave = false
        BBMERGE (
            ch_properly_paired_reads,
            interleave
        )

        // putting together merged and unmerged reads
        ch_processed_paired_reads = BBMERGE.out.merged
                                        .join( BBMERGE.out.unmerged )
                                        .map{
                                            meta, merged, unmerged ->
                                                def reads = [ merged, unmerged ]
                                                [ meta, reads.flatten() ]
                                        }

        ch_versions = ch_versions.mix ( BBMERGE.out.versions )
    }

    // ------------------------------------------------------------------------------------
    // FASTQ TO FASTA
    // ------------------------------------------------------------------------------------

    // putting together single reads and processed paired reads
    ch_reads = ch_branched_reads.single
                    .mix ( ch_processed_paired_reads )

    SEQKIT_FQ2FA ( ch_reads )

    // adding read fasta length to meta
    ch_fasta = SEQKIT_FQ2FA.out.fasta
                    .map {
                        meta, read_fasta_sum_len, file ->
                            [ meta + [ read_fasta_sum_len: read_fasta_sum_len.toLong() ], file ]
                    }

    emit:
    merged_reads_fasta              = ch_fasta
    versions                        = ch_versions

}
