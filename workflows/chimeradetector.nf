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
include { GET_CHIMERAS                                                              } from '../subworkflows/local/get_chimeras'
include { MULTIQC_WORKFLOW                                                          } from '../subworkflows/local/multiqc'

include { NCBI_ASSEMBLY_STATS                                                       } from '../modules/local/ncbi_assembly_stats'
include { SEQKIT_STATS                                                              } from '../modules/nf-core/seqkit/stats'



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
                        "${meta.family}_${meta.taxid}_${meta.sra_id}_${meta.id}",
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
                        "${meta.family}_${meta.taxid}_${meta.sra_id}_${meta.id}",
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
                            [ meta, "${meta.family}_${meta.taxid}_${meta.sra_id}" ]
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
                            .map {
                                meta, files ->
                                    flattened_files = files.flatten()
                                    println "flattened_files $flattened_files"
                                    if ( flattened_files )
                                        [ meta, flattened_files ]
                            }
                            .filter {
                                meta, files -> files.size() == 1 && files[0].toString().contains('to_process')
                            }


    return ch_not_started
            .mix ( ch_started_not_done )
            .map { meta, files -> [ meta, meta.sra_id ] }


}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CHIMERADETECTOR {

    take:
    ch_families
    ch_fastq

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    def target_fasta_file = file( params.target_db, checkExists: true )
    ch_target_db = Channel.of([
        [ id: target_fasta_file.baseName ],
        target_fasta_file
    ])

    // ------------------------------------------------------------------------------------
    // FETCHING STATISTICS FOR ALL FAMILIES INCLUDED IN THE ANALYSIS
    // ------------------------------------------------------------------------------------

    // mix together all families
    ch_all_families = ch_families.mix ( ch_fastq.map { meta, files -> meta.family } )

    NCBI_ASSEMBLY_STATS (
        ch_all_families,
        params.ncbi_api_key ?: []
    )

    // associating back statistics to their respective family channels
    ch_families = ch_families
        .join( NCBI_ASSEMBLY_STATS.out.mean_lengths )
        .map {
            family, mean_assembly_length ->
                def meta = [ family: family, mean_assembly_length: mean_assembly_length ]
                [ meta, family ]
        }

    ch_fastq = ch_fastq
        .combine( NCBI_ASSEMBLY_STATS.out.mean_lengths )
        .filter { meta, files, family, mean_assembly_length -> meta.family == family }
        .map { meta, files, family, mean_assembly_length ->
            def new_meta = meta + [ mean_assembly_length: mean_assembly_length ]
            [ new_meta, files ]
        }

    // ------------------------------------------------------------------------------------
    // GETTING LIST OF SRA IDS
    // ------------------------------------------------------------------------------------

    FETCH_SRA_IDS (
        ch_families,
        params.ncbi_api_key ?: []
    )
    FETCH_SRA_IDS.out.sra_ids.set { ch_sra_ids }
    FETCH_SRA_IDS.out.taxids.set { ch_species_taxids }

    // ------------------------------------------------------------------------------------
    // RESTRICTING TO SPECIFIC SRA IDS / EXCLUDING SPECIFIC SRA IDS
    // ------------------------------------------------------------------------------------

    if ( params.restrict_to_srrs ) {
        restricted_srrs = params.restrict_to_srrs.strip().tokenize(',')
        ch_sra_ids = ch_sra_ids.filter { meta, sra_id -> restricted_srrs.contains(sra_id) }
    }

    if ( params.exclude_srrs ) {
        excluded_srrs = params.exclude_srrs.strip().tokenize(',')
        ch_sra_ids = ch_sra_ids.filter { meta, sra_id -> !excluded_srrs.contains(sra_id) }
    }

    // ------------------------------------------------------------------------------------
    // FILTERING OUT SRA IDS FOR WHICH WE ALREADY HAVE RESULTS
    // TO AVOID PROCESSING THEM AGAINST
    // (IN CASE THE PIPELINE STOPPED AND HAD TO BE RESTARTED)
    // ------------------------------------------------------------------------------------

    ch_sra_ids.set { ch_not_processed_sra_ids }
    if ( !params.ignore_sra_registry ) {
        ch_not_processed_sra_ids = getSraIdsNotProcessed ( ch_sra_ids )
    }

    // ------------------------------------------------------------------------------------
    // DOWNLOAD ALL SRA DATA
    // ------------------------------------------------------------------------------------

    DOWNLOAD_SRA ( ch_not_processed_sra_ids )
    DOWNLOAD_SRA.out.reads.set { ch_sra_reads }

    addToSraRegistry ( ch_sra_reads )

    // mix with custom Fastq files
    ch_sra_reads
        .mix ( ch_fastq )
        .set { ch_reads }

    // ---------------------------------------------------------------
    // COMPUTING STATISTICS FOR EACH SRR / CUSTOM FASTQS
    // ---------------------------------------------------------------

    SEQKIT_STATS ( ch_reads )

    // parsing stats and turning it into a channel readily formated for the MultiQC generalstat table
    SEQKIT_STATS.out.stats
        .splitCsv( header: true, sep: '\t' )
        .map {
            meta, stats_map ->
                ["file", "format", "type"].each { stats_map.remove(it) } // removing fields
                [ id: meta.id ] + stats_map
        }
        .set { ch_fastq_stats }

    // ------------------------------------------------------------------------------------
    // DOWNLOAD GENOME ASSEMBLY FROM NCBI IF AVAILABLE OR MAKE FROM SCRATCH OTHERWISE
    // ------------------------------------------------------------------------------------

    GET_GENOMES (
        ch_reads,
        params.ncbi_api_key ?: []
    )
    GET_GENOMES.out.assemblies.set { ch_assemblies }

    // ------------------------------------------------------------------------------------
    // COMBINE PAIRED READS (IF NECESSARY) AND CONVERT FASTQ TO FASTA
    // ------------------------------------------------------------------------------------

    POST_PROCESS_SRA ( ch_reads )
    POST_PROCESS_SRA.out.merged_reads_fasta.set { ch_reads_fasta }

    // ------------------------------------------------------------------------------------
    // BLAST AGAINST TARGET
    // ------------------------------------------------------------------------------------

    BLAST_AGAINST_TARGET (
        ch_reads_fasta,
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
    // FIND CHIMERAS AND COMPUTE COVERAGE
    // ------------------------------------------------------------------------------------

    GET_CHIMERAS (
        BLAST_AGAINST_TARGET.out.hits,
        BLAST_AGAINST_GENOMES.out.hits
    )
    GET_CHIMERAS.out.chimeras_csv.set { ch_chimeras_csv }

    addDoneToSraRegistry ( ch_chimeras_csv )

    // ------------------------------------------------------------------------------------
    // MULTIQC
    // ------------------------------------------------------------------------------------

    ch_versions
        .mix( DOWNLOAD_SRA.out.versions )
        .mix( SEQKIT_STATS.out.versions )
        .mix( GET_GENOMES.out.versions )
        .mix( POST_PROCESS_SRA.out.versions )
        .mix( BLAST_AGAINST_TARGET.out.versions )
        .mix( BLAST_AGAINST_GENOMES.out.versions )
        .set { ch_versions }


    MULTIQC_WORKFLOW (
        ch_chimeras_csv,
        ch_fastq_stats,
        ch_reads_fasta,
        ch_species_taxids,
        ch_versions
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
