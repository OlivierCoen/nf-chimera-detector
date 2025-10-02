include { BLAST_MAKEBLASTDB as MAKEBLASTDB             } from '../../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN as BLASTN                       } from '../../../modules/local/blast/blastn'


workflow BLAST_AGAINST_GENOMES {

    take:
    ch_sra_reads
    ch_assemblies

    main:

    ch_versions = Channel.empty()

    // ------------------------------------------------------------------------------------
    // BUILD BLAST DB FOR EACH GENOME
    // ------------------------------------------------------------------------------------

    MAKEBLASTDB ( ch_assemblies )

    // ------------------------------------------------------------------------------------
    // BLASTN OF READ FASTA FILES AGAINST THEIR RESPECTIVE GENOME DB
    // ------------------------------------------------------------------------------------

    // associating reads to their respective assembly
    // join by id (SRR ID)
    ch_sra_reads
        .combine( MAKEBLASTDB.out.db )
        .filter { meta, read_fasta, meta_db, db -> meta.id == meta_db.id }
        .map { meta, read_fasta, meta_db, db -> [ meta, read_fasta, db ] }
        .set { blastn_input }

    BLASTN ( blastn_input )


    ch_versions
        .mix ( MAKEBLASTDB.out.versions )
        .set { ch_versions }

    emit:
    hits                            = BLASTN.out.txt
    versions                        = ch_versions

}

