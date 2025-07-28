include { SEQTK_SUBSEQ                                 } from '../../../modules/local/seqtk/subseq'

process EXTRACT_SEQ_IDS {
    input:
    tuple val(meta), path(tsv_file)

    output:
    tuple val(meta), path("${meta.id}_hits.txt"), emit: ids

    script:
    """
    cut -f1 ${tsv_file} > ${meta.id}_hits.txt
    """
}


workflow GET_READS_WITH_HITS {

    take:
    ch_sra_reads
    ch_blast_target_hits

    main:

    ch_versions = Channel.empty()

    EXTRACT_SEQ_IDS ( ch_blast_target_hits )

    ch_sra_reads
        .join ( EXTRACT_SEQ_IDS.out.ids )
        .set { seqtk_subseq_input }

    SEQTK_SUBSEQ ( seqtk_subseq_input )


    emit:
    reads                           = SEQTK_SUBSEQ.out.sequences
    versions                        = ch_versions

}

