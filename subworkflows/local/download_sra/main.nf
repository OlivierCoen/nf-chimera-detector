include { CUSTOM_SRATOOLSNCBISETTINGS } from '../../../modules/nf-core/custom/sratoolsncbisettings'
include { SRATOOLS_PREFETCH           } from '../../../modules/local/sratools/prefetch'
include { SRATOOLS_FASTERQDUMP        } from '../../../modules/local/sratools/fasterqdump'

// ----------------------------------------------------------------------------
// DOWNLOAD FASTQ SEQUENCING READS FROM THE NCBI'S SEQUENCE READ ARCHIVE (SRA).
// ----------------------------------------------------------------------------

workflow DOWNLOAD_SRA {
    take:
    ch_sra_ids   // channel: [ val(meta), val(sra_id) ]

    main:

    ch_versions = Channel.empty()

    // --------------------------------------------------------
    // DETECT EXISTING NCBI USER SETTINGS OR CREATE NEW ONES.
    // --------------------------------------------------------

    CUSTOM_SRATOOLSNCBISETTINGS ( ch_sra_ids.collect() )
    ch_ncbi_settings = CUSTOM_SRATOOLSNCBISETTINGS.out.ncbi_settings

    // ----------------------------------------
    // PREFETCH SEQUENCING READS IN SRA FORMAT.
    // ----------------------------------------

    SRATOOLS_PREFETCH (
        ch_sra_ids,
        ch_ncbi_settings
    )

    SRATOOLS_PREFETCH.out.sra
        .map {
            meta, sra ->
                def new_meta = [ id: sra.name ] + meta
                [ new_meta, sra ]
        }
        .transpose() // when multiple SRRs are downloaded for a specific SRA ID, we split them
        .set { ch_sra }

    // ---------------------------------------------------------------
    // CONVERT THE SRA FORMAT INTO ONE OR MORE COMPRESSED FASTQ FILES.
    // ---------------------------------------------------------------

    SRATOOLS_FASTERQDUMP (
        ch_sra,
        ch_ncbi_settings
    )
    SRATOOLS_FASTERQDUMP.out.reads.set { ch_sra_reads }


    ch_versions
        .mix( CUSTOM_SRATOOLSNCBISETTINGS.out.versions )
        .set { ch_versions }


    emit:
    reads               = ch_sra_reads
    versions            = ch_versions
}
