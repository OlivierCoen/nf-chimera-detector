include { GET_CHILDREN_TAXIDS                                     } from '../../../modules/local/get_children_taxids'
include { GET_SRA_METADATA                                        } from '../../../modules/local/get_sra_metadata'


workflow FETCH_SRA_IDS {

    take:
    ch_families
    ncbi_api_key

    main:

    // ------------------------------------------------------------------------------------
    // GET_CHILDREN_TAXIDS is USELESS NOW
    // WE KEEP IT ONLY FOR THE SYLE (ACTUALLY MULTIQC REPORT)
    // ------------------------------------------------------------------------------------

    GET_CHILDREN_TAXIDS (
        ch_families,
        ncbi_api_key
    )

    ch_species_taxids = GET_CHILDREN_TAXIDS.out.taxid_to_names_files
                            .map { meta, file -> [ meta, file.splitCsv( header: ['taxid', 'taxon_name'] ) ] }
                            .transpose() // explodes each list : we get items like [[family: ..., mean_assembly_length: ...], [taxid: ..., taxon_name: ...]]
                            .map { metas -> metas.collectEntries { it } } // flattens both maps together
                            .map { meta -> [ meta, meta.taxid.strip() ] }

    GET_SRA_METADATA ( ch_families )

    // ------------------------------------------------------------------------------------
    // FOR DEV PURPOSES : RESTRICTING SELECTED SRRS
    // ------------------------------------------------------------------------------------

    ch_sra_ids = GET_SRA_METADATA.out.taxid_sra_id_file
                    .map {
                        meta, taxid, file ->
                            def new_meta = meta + [taxid: taxid]
                            if ( params.max_srrs_per_taxid ) { // in dev, limiting the nb of SRR per taxid
                                [ new_meta, file.splitText( limit: params.max_srrs_per_taxid ) ]
                            } else {
                                [ new_meta, file.splitText() ]
                            }
                    }

    // ------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------

    ch_sra_ids = ch_sra_ids
                    .transpose() // explodes each list
                    .unique() // there may be duplicates
                    .map {
                        meta, sra_id ->
                            def new_meta = meta + [ sra_id: sra_id.strip() ]
                            [ new_meta, sra_id.strip() ]
                    }

    emit:
    sra_ids           = ch_sra_ids
    taxids            = ch_species_taxids
}
