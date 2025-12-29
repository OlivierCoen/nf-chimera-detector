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
    ch_target_hits
    .join( ch_genome_hits )
    .collectFile(
        name: 'test.tsv',
        seed: "id\tfamily\tmean_assembly_length\ttaxid\ttaxon_name\tsra_id\tread_fasta_sum_len\tfile1\tfile2", // header of TSV file
        newLine: true,
        storeDir: "${params.outdir}"
    ){
        meta, f1, f2 -> "${meta.id}\t${meta.family}\t${meta.mean_assembly_length}\ttxid${meta.taxid}\t${meta.taxon_name}\t${meta.sra_id}\t${meta.read_fasta_sum_len}\t${f1}\t${f2}"
    }


    FIND_CHIMERAS (
        ch_target_hits.join( ch_genome_hits )
    )
    FIND_CHIMERAS.out.csv.view { meta, file -> "chimera ${meta}, ${file}" }.set { ch_chimeras_csv }

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
