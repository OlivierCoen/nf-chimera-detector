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
        .map {
            meta, taxid ->
                def new_meta = meta + [ taxid: taxid.strip() ]
                [ new_meta, taxid.strip() ]
        }
        .set { ch_species_taxids }

    GET_SRA_METADATA ( ch_species_taxids )

    GET_SRA_METADATA.out.sra_id_files
        .map {
            meta, file ->
                if ( params.max_srrs_per_taxid ) { // in dev, limiting the nb of SRR per taxid
                    [ meta, file.splitText( limit: params.max_srrs_per_taxid ) ]
                } else {
                    [ meta, file.splitText() ]
                }
        }
        .transpose() // explodes each list
        .map {
            meta, sra_id ->
                def new_meta = meta + [ original_sra_id: sra_id.strip() ]
                [ new_meta, sra_id.strip() ]
        }
        .set { ch_sra_ids }

    emit:
    sra_ids = ch_sra_ids
}

