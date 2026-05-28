include { FIND_CHIMERAS                       } from '../../../modules/local/find_chimeras'
include { MERGE_CHIMERAS                       } from '../../../modules/local/merge_chimeras'
include { GET_CHIMERA_READ_COVERAGE           } from '../../../modules/local/get_chimera_read_coverage'


// ------------------------------------------------------------------------------------
// WORKFLOW
// ------------------------------------------------------------------------------------


workflow GET_CHIMERAS {

    take:
    ch_target_hits
    ch_genome_hits

    main:

    // ------------------------------------------------------------------------------------
    // FIND CHIMERAS FROM BLAST RESULTS
    // ------------------------------------------------------------------------------------

    FIND_CHIMERAS (
        ch_target_hits.combine( ch_genome_hits, by: 0)
    )

    // ------------------------------------------------------------------------------------
    // MERGE CHIMERAS OBTAINED FROM SPLITTED BLAST RESULTS
    // ------------------------------------------------------------------------------------

    MERGE_CHIMERAS( FIND_CHIMERAS.out.csv.groupTuple() )
    ch_chimeras_csv = MERGE_CHIMERAS.out.chimeras

    // ------------------------------------------------------------------------------------
    // COMPUTE COVERAGE OF CHIMERAS ON TARGETS
    // ------------------------------------------------------------------------------------

    // filter out empty chimera files
    ch_nonempty_chimeras_csv = ch_chimeras_csv.filter { meta, file -> file.size() > 0 }

    GET_CHIMERA_READ_COVERAGE (
        ch_target_hits.join( ch_nonempty_chimeras_csv )
    )

    emit:
    chimeras_csv      = ch_chimeras_csv

}
