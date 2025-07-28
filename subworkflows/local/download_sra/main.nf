include { CUSTOM_SRATOOLSNCBISETTINGS } from '../../../modules/nf-core/custom/sratoolsncbisettings'
include { SRATOOLS_PREFETCH           } from '../../../modules/local/sratools/prefetch'
include { SRATOOLS_FASTERQDUMP        } from '../../../modules/local/sratools/fasterqdump'

//
// Download FASTQ sequencing reads from the NCBI's Sequence Read Archive (SRA).
//
workflow DOWNLOAD_SRA {
    take:
    ch_sra_ids   // channel: [ val(meta), val(sra_id) ]

    main:

    ch_versions = Channel.empty()

    //
    // Detect existing NCBI user settings or create new ones.
    //
    CUSTOM_SRATOOLSNCBISETTINGS ( ch_sra_ids.collect() )
    ch_ncbi_settings = CUSTOM_SRATOOLSNCBISETTINGS.out.ncbi_settings
    ch_versions = ch_versions.mix(CUSTOM_SRATOOLSNCBISETTINGS.out.versions)

    //
    // Prefetch sequencing reads in SRA format.
    //
    SRATOOLS_PREFETCH (
        ch_sra_ids,
        ch_ncbi_settings
    )

    SRATOOLS_PREFETCH.out.sra
        .map {
            meta, sra ->
                def new_meta = meta + [ id: sra.name ]
                [ new_meta, sra ]
        }
        .set { ch_sra }

    //
    // Convert the SRA format into one or more compressed FASTQ files.
    //
    SRATOOLS_FASTERQDUMP (
        ch_sra,
        ch_ncbi_settings
    )

    emit:
    reads    = SRATOOLS_FASTERQDUMP.out.reads // channel: [ val(meta), [ reads ] ]
    versions = ch_versions                    // channel: [ versions.yml ]
}
