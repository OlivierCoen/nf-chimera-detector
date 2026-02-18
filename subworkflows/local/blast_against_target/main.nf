include { BLAST_MAKEBLASTDB as MAKEBLASTDB             } from '../../../modules/local/blast/makeblastdb'
include { BLAST_BLASTN as BLASTN                       } from '../../../modules/local/blast/blastn'
include { SEQTK_SUBSEQ                                 } from '../../../modules/local/seqtk/subseq'
include { SEQKIT_SPLIT2                                } from '../../../modules/local/seqkit/split2/'

// ------------------------------------------------------------------------------------
// PROCESSES
// ------------------------------------------------------------------------------------

process MERGE_HITS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(hit_files, stageAs: "**/*")

    output:
    tuple val(meta), path("${meta.id}.all_hits.txt"), emit: hits

    script:
    """
    cat ${hit_files} > ${meta.id}.all_hits.txt
    """
}


process EXTRACT_SEQ_IDS {
    tag "${meta.id}"
    input:
    tuple val(meta), path(tsv_file)

    output:
    tuple val(meta), path("${meta.id}_hit_ids.txt"), emit: ids

    script:
    """
    cut -f1 ${tsv_file} > ${meta.id}_hit_ids.txt
    """
}


// ------------------------------------------------------------------------------------
// WORKFLOW
// ------------------------------------------------------------------------------------

workflow BLAST_AGAINST_TARGET {

    take:
    ch_reads_fasta
    ch_target_db

    main:

    ch_versions = channel.empty()

    // ------------------------------------------------------------------------------------
    // MAKE DB
    // ------------------------------------------------------------------------------------

    MAKEBLASTDB ( ch_target_db )

    // ------------------------------------------------------------------------------------
    // COMPUTING THE NB OF CHUNKS FOR EACH READ FASTA FILE
    // ------------------------------------------------------------------------------------

    ch_split_input = ch_reads_fasta
                        .map { // computing the nb of chunks necessary
                            meta, fasta ->
                                def ratio = meta.read_fasta_sum_len.toLong() / params.read_fasta_chunk_max_size
                                def nb_chunks = ratio.toInteger() + 1
                                [ meta, fasta, nb_chunks ]
                        }

    // ------------------------------------------------------------------------------------
    // SPLITTING QUERY FASTA FILES IN CHUNKS TO AVOID ISSUES
    // ------------------------------------------------------------------------------------

    SEQKIT_SPLIT2 ( ch_split_input )

    // ------------------------------------------------------------------------------------
    // BLASTN
    // ------------------------------------------------------------------------------------

     // making all combinations of splitted reads + target db
    blastn_input = SEQKIT_SPLIT2.out.reads // list of splitted reads
                    .transpose() // separate splitted reads, each with their own meta
                    .combine( MAKEBLASTDB.out.db )
                    .map { meta, reads, meta2, db ->  [ meta, reads, db ] }

    BLASTN ( blastn_input )


    // ------------------------------------------------------------------------------------
    // REGROUP HITS BY META
    // ------------------------------------------------------------------------------------

    MERGE_HITS( BLASTN.out.txt.groupTuple() )
    ch_hits = MERGE_HITS.out.hits

    // ------------------------------------------------------------------------------------
    // EXTRACT SEQ IDS CORRESPONDING TO HITS
    // ------------------------------------------------------------------------------------

    EXTRACT_SEQ_IDS ( ch_hits )

    seqtk_subseq_input = ch_reads_fasta
                            .join ( EXTRACT_SEQ_IDS.out.ids )


    // ------------------------------------------------------------------------------------
    // GET READS CORRESPONDING TO THESE SEQ IDS
    // ------------------------------------------------------------------------------------

    SEQTK_SUBSEQ ( seqtk_subseq_input )

    ch_versions = ch_versions
                    .mix ( MAKEBLASTDB.out.versions )

    emit:
    hits                            = ch_hits
    hit_sequences                   = SEQTK_SUBSEQ.out.sequences
    versions                        = ch_versions

}
