include { GET_BEST_NCBI_ASSEMBLY                     } from '../../../modules/local/get_best_ncbi_assembly'
include { BLAST_BLASTN                          } from '../../../modules/local/blast/blastn'


workflow ASSEMBLY {

    take:
    ch_sra_reads

    main:

    ch_versions = Channel.empty()

    ch_sra_reads
        .map { meta, reads -> [ meta, meta.taxid ]}
        .set { ch_taxids }

    GET_BEST_NCBI_ASSEMBLY ( ch_taxids )

    ch_sra_reads
        .join ( GET_BEST_NCBI_ASSEMBLY.out.accession )
        .branch {
            meta, reads, accession ->
                to_download: accession != 'NONE'
                to_assemble: accession == 'NONE'
        }
        .set { ch_branched_sra_reads }



    emit:
    hits                            = BLAST_BLASTN.out.txt
    versions                        = ch_versions

}

