include { FIND_CHIMERAS                       } from '../../../modules/local/find_chimeras'
include { GET_CHIMERA_READ_COVERAGE           } from '../../../modules/local/get_chimera_read_coverage'


workflow GET_CHIMERAS {

    take:
    ch_target_hits
    ch_genome_hits

    main:

    // ------------------------------------------------------------------------------------
    // FIND CHIMERAS FROM BLAST RESULTS
    // ------------------------------------------------------------------------------------

    FIND_CHIMERAS (
        ch_target_hits.join( ch_genome_hits )
    )
    FIND_CHIMERAS.out.csv.set { ch_chimeras_csv }

    // ------------------------------------------------------------------------------------
    // COMPUTE COVERAGE OF CHIMERAS ON TARGETS
    // ------------------------------------------------------------------------------------

    // filter out empty chimera files
    ch_chimeras_csv
        .filter { meta, file -> file.size() > 0 }
        .set { ch_nonempty_chimeras_csv }

    GET_CHIMERA_READ_COVERAGE (
        ch_target_hits.join( ch_nonempty_chimeras_csv )
    )

    emit:
    chimeras_csv      = ch_chimeras_csv

}
