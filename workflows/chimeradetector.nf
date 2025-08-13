/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FETCH_SRA_IDS                                                             } from '../subworkflows/local/fetch_sra_ids'
include { DOWNLOAD_SRA                                                              } from '../subworkflows/local/download_sra'
include { GET_GENOMES                                                               } from '../subworkflows/local/get_genomes'
include { POST_PROCESS_SRA                                                          } from '../subworkflows/local/post_process_sra'
include { BLAST_AGAINST_TARGET                                                      } from '../subworkflows/local/blast_against_target'
include { BLAST_AGAINST_GENOMES                                                     } from '../subworkflows/local/blast_against_genomes'

include { NCBI_ASSEMBLY_STATS                                                       } from '../modules/local/ncbi_assembly_stats'
include { FIND_CHIMERAS                                                             } from '../modules/local/find_chimeras'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// STORING SRA IDS TO PROCESS
def addToSraRegistry ( ch_sra_reads) {
    ch_sra_reads
        .collectFile (
            {
                meta, reads ->
                    [
                        "${meta.family}_${meta.taxid}_${meta.original_sra_id}_${meta.id}",
                        reads instanceof List ? reads.collect { it.name }.join("\n") : reads.name
                    ]
            },
            storeDir: "${params.outdir}/${params.sra_registry}/to_process/"
        )
}

// STORING SRA IDS DONE
def addDoneToSraRegistry ( ch_csv ) {
    ch_csv
        .collectFile (
            {
                meta, csv ->
                    [
                        "${meta.family}_${meta.taxid}_${meta.original_sra_id}_${meta.id}",
                        csv.name
                    ]
            },
            storeDir: "${params.outdir}/${params.sra_registry}/done/"
        )
}

def parseSraFolder ( ch_prefix, folder ) {
    return ch_prefix
                .map {
                    meta, prefix ->
                        def matchingFiles = folder
                                                .listFiles()
                                                .findAll { it.name.startsWith(prefix) }
                        [ meta, matchingFiles ]
                }
                .groupTuple()
}


def groupFilesBySRR ( fileList ) {
    // group by SRR ID
    return fileList
                .groupBy { it.name.split('_').last() }
                .collect { prefix, files -> [ files ] }
}

// GET LIST OF SRA IDS THAT ALL NOT ALL DONE
// A SRA ID IS CONSIDERED DONE WHEN ALL ITS UNDERLYING SRRS ARE MARKED AS DONE
def getSraIdsNotProcessed ( ch_sra_ids ) {

    ch_prefix = ch_sra_ids
                    .map {
                        meta, _ ->
                            [ meta, "${meta.family}_${meta.taxid}_${meta.original_sra_id}" ]
                    }

    def to_process_folder = file ( "${params.outdir}/${params.sra_registry}/to_process/" )
    def done_folder = file ( "${params.outdir}/${params.sra_registry}/done/" )

    ch_to_process = parseSraFolder ( ch_prefix, to_process_folder )
    ch_done = parseSraFolder ( ch_prefix, done_folder )

    ch_grouped = ch_to_process
                    .mix( ch_done )
                    .groupTuple()
                    .map { meta, files -> [ meta, files.flatten() ] }

    ch_started = ch_grouped.filter { meta, files -> files.size() > 0 }
    ch_not_started = ch_grouped.filter { meta, files -> files.size() == 0 }

    ch_started_not_done = ch_started
                            .map { meta, files -> [ meta, groupFilesBySRR ( files ) ] }
                            .transpose()
                            .map { meta, files -> [ meta, files.flatten() ] }
                            .filter {
                                meta, files -> files.size() == 1 && files[0].toString().contains('to_process')
                            }


    return ch_not_started
            .mix ( ch_started_not_done )
            .map { meta, files -> [ meta, meta.original_sra_id ] }


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHIMERADETECTOR {

    take:
    ch_families

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    def target_fasta_file = file( params.target_db, checkExists: true )
    ch_target_db = Channel.of([
        [ id: target_fasta_file.baseName ],
        target_fasta_file
    ])

    NCBI_ASSEMBLY_STATS ( ch_families )

    NCBI_ASSEMBLY_STATS.out.mean_lengths
        .map {
            family, mean_assembly_length ->
                def meta = [ family: family, mean_assembly_length: mean_assembly_length ]
                [ meta, family ]
        }
        .set { ch_families }

    // ------------------------------------------------------------------------------------
    // GETTING LIST OF SRA IDS
    // ------------------------------------------------------------------------------------

    FETCH_SRA_IDS ( ch_families )
    FETCH_SRA_IDS.out.sra_ids.set { ch_sra_ids }

    // For dev purposes: test download of a specific SRA accession
    if ( params.only_download_sra ) {

        unique_sra_accession = params.only_download_sra.strip()
        ch_sra_ids = ch_sra_ids.filter { meta, sra_id -> sra_id == unique_sra_accession }

    }

    // ------------------------------------------------------------------------------------
    // FILTERING OUT SRA IDS FOR WHICH WE ALREADY HAVE RESULTS
    // TO AVOID PROCESSING THEM AGAINST
    // (IN CASE THE PIPELINE STOPPED AND HAD TO BE RESTARTED)
    // ------------------------------------------------------------------------------------

    ch_not_processed_sra_ids = getSraIdsNotProcessed ( ch_sra_ids )

    // ------------------------------------------------------------------------------------
    // DOWNLOAD ALL SRA DATA
    // ------------------------------------------------------------------------------------

    DOWNLOAD_SRA ( ch_not_processed_sra_ids )
    DOWNLOAD_SRA.out.reads.set { ch_sra_reads }

    addToSraRegistry ( ch_sra_reads )

    // ------------------------------------------------------------------------------------
    // DOWNLOAD GENOME ASSEMBLY FROM NCBI IF AVAILABLE OR MAKE FROM SCRATCH OTHERWISE
    // ------------------------------------------------------------------------------------

    GET_GENOMES  ( ch_sra_reads )
    GET_GENOMES .out.assemblies.set { ch_assemblies }

    // ------------------------------------------------------------------------------------
    // COMBINE PAIRED READS (IF NECESSARY) AND CONVERT FASTQ TO FASTA
    // ------------------------------------------------------------------------------------

    POST_PROCESS_SRA ( ch_sra_reads )
    POST_PROCESS_SRA.out.single_reads.set { ch_sra_reads }

    // ------------------------------------------------------------------------------------
    // BLAST AGAINST TARGET
    // ------------------------------------------------------------------------------------

    BLAST_AGAINST_TARGET (
        ch_sra_reads,
        ch_target_db
    )

    // ------------------------------------------------------------------------------------
    // BLAST AGAINST GENOMES (DOWNLOADED OR ASSEMBLED)
    // ------------------------------------------------------------------------------------

    BLAST_AGAINST_GENOMES (
        BLAST_AGAINST_TARGET.out.hit_sequences,
        ch_assemblies
    )

    // ------------------------------------------------------------------------------------
    // FIND CHIMERAS
    // ------------------------------------------------------------------------------------

    BLAST_AGAINST_TARGET.out.hits
        .join( BLAST_AGAINST_GENOMES.out.hits )
        .set { find_chimeras_input }

    FIND_CHIMERAS ( find_chimeras_input )

    addDoneToSraRegistry ( FIND_CHIMERAS.out.csv )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
