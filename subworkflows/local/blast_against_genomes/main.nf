include { BLAST_MAKEBLASTDB as MAKEBLASTDB             } from '../../../modules/local/blast/makeblastdb'
include { BLAST_BLASTN as BLASTN                       } from '../../../modules/local/blast/blastn'
include { SEQTK_SUBSEQ                                 } from '../../../modules/local/seqtk/subseq'


workflow BLAST_AGAINST_GENOMES {

    take:
    ch_reads_fasta
    ch_hit_ids
    ch_assemblies

    main:

    ch_versions = channel.empty()

    // ------------------------------------------------------------------------------------
    // BUILD BLAST DB FOR EACH GENOME
    // ------------------------------------------------------------------------------------

    MAKEBLASTDB ( ch_assemblies )

    // ------------------------------------------------------------------------------------
    // GET READS CORRESPONDING TO THESE SEQ IDS
    // ------------------------------------------------------------------------------------

    SEQTK_SUBSEQ (
        ch_reads_fasta.join( ch_hit_ids )
    )

    // ------------------------------------------------------------------------------------
    // BLASTN OF READ FASTA FILES AGAINST THEIR RESPECTIVE GENOME DB
    // ------------------------------------------------------------------------------------

    // associating reads to their respective assembly
    // join by id (SRR ID)
    blastn_input = SEQTK_SUBSEQ.out.sequences
                    .combine( MAKEBLASTDB.out.db )
                    .filter { meta, read_fasta, meta_db, db -> meta.id == meta_db.id }
                    .map { meta, read_fasta, meta_db, db -> [ meta, read_fasta, db ] }

    BLASTN ( blastn_input )


    ch_versions = ch_versions
        .mix ( MAKEBLASTDB.out.versions )

    emit:
    hits                            = BLASTN.out.txt
    versions                        = ch_versions

}
