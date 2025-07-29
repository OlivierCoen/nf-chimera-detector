include { BLAST_MAKEBLASTDB as MAKEBLASTDB             } from '../../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN as BLASTN                       } from '../../../modules/local/blast/blastn'
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


workflow BLAST_AGAINST_TARGET {

    take:
    ch_sra_reads
    ch_target_db

    main:

    ch_versions = Channel.empty()

    // ------------------------------------------------------------------------------------
    // MAKE DB
    // ------------------------------------------------------------------------------------

    MAKEBLASTDB ( ch_target_db )

    // ------------------------------------------------------------------------------------
    // BLASTN
    // ------------------------------------------------------------------------------------

     // making all combinations of reads + target db
    ch_sra_reads
        .combine( MAKEBLASTDB.out.db )
        .map { meta, reads, meta2, db ->  [ meta, reads, db ] }
        .set { blastn_input }

    BLASTN ( blastn_input )

    BLASTN.out.txt.set { ch_hits }

    // ------------------------------------------------------------------------------------
    // EXTRACT SEQ IDS CORRESPONDING TO HITS
    // ------------------------------------------------------------------------------------

    EXTRACT_SEQ_IDS ( ch_hits )

    ch_sra_reads
        .join ( EXTRACT_SEQ_IDS.out.ids )
        .set { seqtk_subseq_input }

    // ------------------------------------------------------------------------------------
    // GET READS CORRESPONDING TO THESE SEQ IDS
    // ------------------------------------------------------------------------------------

    SEQTK_SUBSEQ ( seqtk_subseq_input )

    ch_versions
        .mix ( MAKEBLASTDB.out.versions )
        .set { ch_versions }

    emit:
    hits                            = ch_hits
    hit_sequences                   = SEQTK_SUBSEQ.out.sequences
    versions                        = ch_versions

}

