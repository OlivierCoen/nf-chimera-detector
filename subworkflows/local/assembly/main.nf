include { GET_BEST_NCBI_ASSEMBLY                     } from '../../../modules/local/get_best_ncbi_assembly'
include { DOWNLOAD_NCBI_ASSEMBLY                     } from '../../../modules/local/download_ncbi_assembly'
include { MEGAHIT                                    } from '../../../modules/nf-core/megahit'


workflow ASSEMBLY {

    take:
    ch_sra_reads

    main:

    ch_versions = Channel.empty()

    ch_sra_reads
        .map { meta, reads -> [ meta, meta.taxid ]}
        .set { ch_taxids }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FETCH ACCESSIONS OF AVAILABLE ASSEMBLIES FROM NCBI
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GET_BEST_NCBI_ASSEMBLY ( ch_taxids )

    ch_sra_reads
        .join ( GET_BEST_NCBI_ASSEMBLY.out.accession )
        .branch {
            meta, reads, accession ->
                to_download: accession != 'NONE'
                to_assemble: accession == 'NONE'
        }
        .set { ch_branched_sra_reads }

    ch_branched_sra_reads.to_download
        .map { meta, reads, accession -> [ meta, accession ] }
        .unique()
        .set { ch_accessions_to_download }

    ch_branched_sra_reads.to_assemble
        .map { meta, reads, accession -> [ meta, reads ] }
        .set { ch_sra_reads_to_assemble }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DOWNLOAD AVAILABLE ASSEMBLIES
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DOWNLOAD_NCBI_ASSEMBLY ( ch_accessions_to_download )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ASSEMBLE GENOMES WHENEVER NO ASSEMBLY IS AVAILABLE ON NCBI
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_sra_reads_to_assemble
        .view { v -> "to assemble" + v }


    MEGAHIT ( ch_sra_reads_to_assemble )

    DOWNLOAD_NCBI_ASSEMBLY.out.assemblies
        .mix ( MEGAHIT.out.contigs )
        .set { ch_assemblies }

    ch_versions = ch_versions
                    .mix ( MEGAHIT.out.versions )

    emit:
    assemblies                      = ch_assemblies
    versions                        = ch_versions

}

