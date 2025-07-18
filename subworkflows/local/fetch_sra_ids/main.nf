include { GET_SPECIES_TAXIDS                                      } from '../../../modules/local/get_species_taxids'
include { GET_SRA_METADATA                                        } from '../../../modules/local/get_sra_metadata'


workflow FETCH_SRA_IDS {

    take:
    ch_families

    main:

    ch_versions = Channel.empty()

    GET_SPECIES_TAXIDS ( ch_families )
    GET_SPECIES_TAXIDS.out.taxid_files
    .map
        .splitText()
        .map { taxid -> taxid.strip() }
        .set { ch_species_taxids }

    GET_SRA_METADATA ( ch_species_taxids )
    GET_SRA_METADATA.out.sra_id_files
        .splitText()
        .map { sra_id -> sra_id.strip() }
        .set { ch_sra_ids }

    emit:
    sra_ids = ch_sra_ids
}

