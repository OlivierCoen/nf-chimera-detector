include { GET_BEST_NCBI_ASSEMBLY                     } from '../../../modules/local/get_best_ncbi_assembly'
include { DOWNLOAD_NCBI_ASSEMBLY                     } from '../../../modules/local/download_ncbi_assembly'
include { MEGAHIT                                    } from '../../../modules/local/megahit'

include { SRA_READS_PREPARATION                      } from '../sra_reads_preparation'


workflow GET_GENOMES {

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

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DOWNLOAD AVAILABLE ASSEMBLIES
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_branched_sra_reads.to_download
        .map { meta, reads, accession -> [ [ family: meta.family, taxid: meta.taxid ], accession ] }
        .unique()
        .set { accessions_to_download }

    DOWNLOAD_NCBI_ASSEMBLY ( accessions_to_download )

    // associating downloaded assemblies to all possible corresponding SRRs
    // there can be multiple SRRs for one single taxid (and hence one single downloaded assembly)
    // in this case, we create a channel with one assembly per SRR
    // we link the assemblies with their SRRs through meta
    DOWNLOAD_NCBI_ASSEMBLY.out.assemblies
        .map { meta, assembly -> [ meta, assembly, meta.taxid ] }
        .set { ch_downloaded_assemblies }

    ch_branched_sra_reads.to_download
        .map { meta, reads, accession -> [ meta, reads, meta.taxid ] }
        .combine( ch_downloaded_assemblies, by: 2 ) // combining by taxid
        .map { taxid, meta_reads, reads, meta_assembly, assembly -> [ meta_reads + meta_assembly, assembly ] } // combining meta maps
        .set { ch_downloaded_assemblies }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ASSEMBLE GENOMES WHENEVER NO ASSEMBLY IS AVAILABLE ON NCBI
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // SUBSAMPLE AND CLEAN READS
    ch_branched_sra_reads.to_assemble
        .map { meta, reads, accession -> [ meta, reads ] }
        .set { ch_sra_reads_to_prepare }

    SRA_READS_PREPARATION ( ch_sra_reads_to_prepare )

    // arranging channel : some are single reads and some other are paired reads
    SRA_READS_PREPARATION.out.prepared_sra_reads
        .map {
            meta, reads ->
                if ( reads instanceof Path ) {
                    [ meta, reads, [] ]
                } else { // List
                    def ( reads_1, reads_2 ) = reads
                    [ meta, reads_1, reads_2 ]
                }
        }
        .set { ch_sra_reads_to_assemble }

    // ASSEMBLE
    MEGAHIT ( ch_sra_reads_to_assemble )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GATHERING ALL ASSEMBLIES TOGETHER
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

     ch_downloaded_assemblies
        .mix ( MEGAHIT.out.contigs )
        .set { ch_assemblies }

    ch_versions = ch_versions
                    .mix ( SRA_READS_PREPARATION.out.versions )

    emit:
    assemblies                      = ch_assemblies
    versions                        = ch_versions

}

