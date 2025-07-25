include { SEQKIT_STATS                                              } from '../../../modules/nf-core/seqkit/stats'
include { SEQTK_SAMPLE                                              } from '../../../modules/local/seqtk/sample/main'
include { FASTP                                                     } from '../../../modules/nf-core/fastp'


workflow SRA_READS_PREPARATION {

    take:
    ch_sra_reads

    main:

    ch_versions = Channel.empty()

    // ---------------------------------------------------------------------
    // Statistics from read files (in particular, total nb of bases)
    // ---------------------------------------------------------------------

    SEQKIT_STATS ( ch_sra_reads )

    // ---------------------------------------------------------------------
    // Subsampling reads in read files
    // ---------------------------------------------------------------------

    SEQKIT_STATS.out.stats
        .splitCsv( header: true, sep: '\t', limit: 1 )
        .map { meta, row -> [ meta, row.sum_len as Float ] }
        .join ( ch_sra_reads )
        .map {
            meta, nb_reads, reads ->
                def observed_coverage = nb_reads / meta.mean_assembly_length.toFloat()
                def fraction_to_keep = params.read_coverage / observed_coverage
                [ meta, reads, fraction_to_keep ]
        }
        .set { seqtk_sample_input }

    SEQTK_SAMPLE ( seqtk_sample_input )

    // ---------------------------------------------------------------------
    // Trimming / Filtering
    // ---------------------------------------------------------------------

    FASTP (
        SEQTK_SAMPLE.out.reads,
        [], false, false, true
    )

    ch_versions = ch_versions
                     .mix ( FASTP.out.versions )
                     .mix ( SEQKIT_STATS.out.versions )

    emit:
    prepared_sra_reads              = FASTP.out.reads
    versions                        = ch_versions
}
