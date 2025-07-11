/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC_WORKFLOW                                        } from '../subworkflows/local/multiqc'

include { GET_SPECIES_TAXIDS                                      } from '../modules/local/get_species_taxids'
include { GET_SRRS                                                } from '../modules/local/get_srrs'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VIRUSLTEFINDER {

    take:
    ch_families

    main:

    ch_versions = Channel.empty()

    GET_SPECIES_TAXIDS ( ch_families )
    GET_SPECIES_TAXIDS.out.taxids
        .splitText()
        .map { taxid -> taxid.strip() }
        .set { ch_species_taxids }

    GET_SRRS ( ch_species_taxids )
    GET_SRRS.out.srrs.view { srr -> "srr " + srr}

    // ------------------------------------------------------------------------------------
    // MULTIQC
    // ------------------------------------------------------------------------------------

    MULTIQC_WORKFLOW ( ch_versions )

    emit:
    multiqc_report = MULTIQC_WORKFLOW.out.multiqc_report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
