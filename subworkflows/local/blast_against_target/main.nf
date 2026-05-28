include { BLAST_MAKEBLASTDB as MAKEBLASTDB             } from '../../../modules/local/blast/makeblastdb'
include { BLAST_BLASTN as BLASTN                       } from '../../../modules/local/blast/blastn'

include { SEQKIT_SPLIT2                                } from '../../../modules/local/seqkit/split2/'

// ------------------------------------------------------------------------------------
// PROCESSES
// ------------------------------------------------------------------------------------


process EXTRACT_SEQ_IDS {
    tag "${meta.id}"
    input:
    tuple val(meta), path(hit_files, stageAs: "**/*")

    output:
    tuple val(meta), path("${meta.id}_hit_ids.txt"), emit: ids

    script:
    """
    for hit_file in ${hit_files}; do
        echo "Extracting seq ids from \$hit_file"
        cut -f1 \$hit_file | sort | uniq >> tmp.txt
    done
    echo "Merging seq ids"
    cat tmp.txt | sort | uniq > ${meta.id}_hit_ids.txt
    rm tmp.txt
    """
}


// ------------------------------------------------------------------------------------
// WORKFLOW
// ------------------------------------------------------------------------------------

workflow BLAST_AGAINST_TARGET {

    take:
    ch_reads_fasta
    ch_target_db
    read_fasta_chunk_max_size

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
                                def ratio = meta.read_fasta_sum_len.toLong() / read_fasta_chunk_max_size
                                def nb_chunks = ratio.toInteger() + 1
                                [ meta, fasta, nb_chunks ]
                        }

    // ------------------------------------------------------------------------------------
    // SPLITTING QUERY FASTA FILES IN CHUNKS TO AVOID ISSUES
    // ------------------------------------------------------------------------------------

    SEQKIT_SPLIT2 ( ch_split_input )
    ch_splitted_read_fasta = SEQKIT_SPLIT2.out.reads

    // ------------------------------------------------------------------------------------
    // BLASTN
    // ------------------------------------------------------------------------------------

     // making all combinations of splitted reads + target db
    blastn_input = ch_splitted_read_fasta // list of splitted reads
                    .transpose() // separate splitted reads, each with their own meta
                    .combine( MAKEBLASTDB.out.db )
                    .map { meta, reads, meta2, db ->  [ meta, reads, db ] }

    BLASTN ( blastn_input )
    ch_hits = BLASTN.out.txt

    // ------------------------------------------------------------------------------------
    // EXTRACT SEQ IDS CORRESPONDING TO HITS
    // ------------------------------------------------------------------------------------

    EXTRACT_SEQ_IDS ( ch_hits.groupTuple() )


    ch_versions = ch_versions
                    .mix ( MAKEBLASTDB.out.versions )

    emit:
    hits                            = ch_hits
    hit_ids                         = EXTRACT_SEQ_IDS.out.ids
    versions                        = ch_versions

}
