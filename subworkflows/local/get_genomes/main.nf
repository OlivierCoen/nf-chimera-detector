include { GET_BEST_NCBI_ASSEMBLY                     } from '../../../modules/local/get_best_ncbi_assembly'
include { DOWNLOAD_NCBI_ASSEMBLY                     } from '../../../modules/local/download_ncbi_assembly'
include { MEGAHIT                                    } from '../../../modules/local/megahit'

include { READ_PREPARATION                           } from '../read_preparation'


workflow GET_GENOMES {

    take:
    ch_sra_reads
    ncbi_api_key

    main:

    ch_versions = channel.empty()

    ch_taxids = ch_sra_reads
                    .map { meta, reads -> [ meta, meta.taxid ]}

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // FETCH ACCESSIONS OF AVAILABLE ASSEMBLIES FROM NCBI
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GET_BEST_NCBI_ASSEMBLY (
        ch_taxids,
        ncbi_api_key
    )

    ch_branched_sra_reads = ch_sra_reads
                                .join ( GET_BEST_NCBI_ASSEMBLY.out.accession )
                                .branch {
                                    meta, reads, accession ->
                                        to_download: accession != 'NONE'
                                        to_assemble: accession == 'NONE'
                                }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // DOWNLOAD AVAILABLE ASSEMBLIES
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    accessions_to_download = ch_branched_sra_reads.to_download
                                .map {
                                    meta, reads, accession -> // remove reads and insert genome_accession in meta map
                                        [ [ family: meta.family, taxid: meta.taxid, genome_accession: accession ], accession ]
                                }
                                .unique()

    DOWNLOAD_NCBI_ASSEMBLY ( accessions_to_download )

    // associating downloaded assemblies to all possible corresponding SRRs
    // there can be multiple SRRs for one single taxid (and hence one single downloaded assembly)
    // in this case, we create a channel with one assembly per SRR
    // we link the assemblies with their SRRs through meta
    ch_downloaded_assemblies = DOWNLOAD_NCBI_ASSEMBLY.out.assemblies
                                .map {
                                    meta, assembly ->
                                        [ meta + [ assembly_name: assembly.baseName ], assembly, meta.taxid ]
                                }

    ch_downloaded_assemblies = ch_branched_sra_reads.to_download
                                    .map { meta, reads, accession -> [ meta, reads, meta.taxid ] }
                                    .combine( ch_downloaded_assemblies, by: 2 ) // combining by taxid
                                    .map { // adding meta maps together
                                        taxid, meta_reads, reads, meta_assembly, assembly ->
                                            [ meta_reads + meta_assembly, assembly ]
                                    }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ASSEMBLE GENOMES WHENEVER NO ASSEMBLY IS AVAILABLE ON NCBI
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // SUBSAMPLE AND CLEAN READS
    ch_reads_to_prepare = ch_branched_sra_reads.to_assemble
                                .map { meta, reads, accession -> [ meta, reads ] }

    READ_PREPARATION ( ch_reads_to_prepare )

    // arranging channel : some are single reads and some other are paired reads
    ch_sra_reads_to_assemble = READ_PREPARATION.out.prepared_sra_reads
                                    .map {
                                        meta, reads ->
                                            if ( reads instanceof Path ) {
                                                [ meta, reads, [] ]
                                            } else { // List
                                                def ( reads_1, reads_2 ) = reads
                                                [ meta, reads_1, reads_2 ]
                                            }
                                    }

    // ASSEMBLE
    MEGAHIT ( ch_sra_reads_to_assemble )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // GATHERING ALL ASSEMBLIES TOGETHER
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

     ch_assemblies = ch_downloaded_assemblies
                        .mix ( MEGAHIT.out.contigs )

    ch_versions = ch_versions
                    .mix ( READ_PREPARATION.out.versions )

    emit:
    assemblies                      = ch_assemblies
    versions                        = ch_versions

}
