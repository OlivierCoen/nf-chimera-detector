include { BLAST_MAKEBLASTDB as MAKEBLASTDB             } from '../../../modules/nf-core/blast/makeblastdb'
include { BLAST_BLASTN as BLASTN                       } from '../../../modules/local/blast/blastn'


workflow BLAST_AGAINST_ASSEMBLIES {

    take:
    ch_sra_reads
    ch_assemblies

    main:

    ch_versions = Channel.empty()

    MAKEBLASTDB ( ch_assemblies )

    // associating reads to their respective assembly
    ch_sra_reads
        .join( MAKEBLASTDB.out.db )
        .set { blastn_input }

    BLASTN ( blastn_input )

    ch_versions
        .mix ( MAKEBLASTDB.out.versions )
        .set { ch_versions }

    emit:
    hits                            = BLASTN.out.txt
    versions                        = ch_versions

}

