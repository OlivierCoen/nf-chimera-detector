include { GET_BEST_NCBI_ASSEMBLY                     } from '../../../modules/local/get_best_ncbi_assembly'
include { DOWNLOAD_NCBI_ASSEMBLY                     } from '../../../modules/local/download_ncbi_assembly'
include { MEGAHIT                                    } from '../../../modules/nf-core/megahit'

include { SRA_READS_PREPARATION                      } from '../sra_reads_preparation'


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
                    [ meta, accession ]
                to_assemble: accession == 'NONE'
                    [ meta, reads ]
        }
        .set { ch_branched_sra_reads }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DOWNLOAD AVAILABLE ASSEMBLIES
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DOWNLOAD_NCBI_ASSEMBLY ( ch_branched_sra_reads.to_download.unique() )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ASSEMBLE GENOMES WHENEVER NO ASSEMBLY IS AVAILABLE ON NCBI
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    SRA_READS_PREPARATION ( ch_sra_reads )

    SRA_READS_PREPARATION.out.prepared_sra_reads
        .map {
            meta, reads ->
                if ( reads instanceof Path ) {
                    [ meta, reads, [] ]
                } else {
                    def ( reads_1, reads_2 ) = reads
                    [ meta, reads_1, reads_2 ]
                }
        }
        .set { ch_sra_reads_to_assemble }

    /*
    MEGAHIT ( ch_sra_reads_to_assemble )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GATHERING ALL ASSEMBLIES TOGETHER
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DOWNLOAD_NCBI_ASSEMBLY.out.assemblies
        .mix ( MEGAHIT.out.contigs )
        .set { ch_assemblies }

    ch_versions = ch_versions
                    .mix ( MEGAHIT.out.versions )
    */
    ch_assemblies = Channel.empty()
    emit:
    assemblies                      = ch_assemblies
    versions                        = ch_versions

}

