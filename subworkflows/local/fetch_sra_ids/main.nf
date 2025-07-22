include { GET_CHILDREN_TAXIDS                                     } from '../../../modules/local/get_children_taxids'
include { GET_SRA_METADATA                                        } from '../../../modules/local/get_sra_metadata'


workflow FETCH_SRA_IDS {

    take:
    ch_families

    main:

    GET_CHILDREN_TAXIDS ( ch_families )
    GET_CHILDREN_TAXIDS.out.taxid_files
        .map { meta, file -> [ meta, file.splitText() ] }
        .transpose() // explodes each list
        .map { meta, taxid -> [ meta, taxid.strip() ] }
        .set { ch_species_taxids }

    GET_SRA_METADATA ( ch_species_taxids )
    GET_SRA_METADATA.out.sra_id_files
        .map { meta, file -> [ meta, file.splitText() ] }
        .transpose() // explodes each list
        .map { meta, sra_id -> [ meta, sra_id.strip() ] }
        .set { ch_sra_ids }

    emit:
    sra_ids = ch_sra_ids
}

