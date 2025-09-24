include { MULTIQC                                                           } from '../../../modules/nf-core/multiqc'
include { PREPARE_DATA_PER_FAMILY as PREPARE_DOWNLOADED_GENOME_SIZE         } from '../../../modules/local/prepare_data_per_family'
include { PREPARE_DATA_PER_FAMILY as PREPARE_ASSEMBLED_GENOME_SIZE          } from '../../../modules/local/prepare_data_per_family'
include { PREPARE_DATA_PER_FAMILY as PREPARE_BLAST_HIT_TARGET               } from '../../../modules/local/prepare_data_per_family'
include { PREPARE_DATA_PER_FAMILY as PREPARE_BLAST_HIT_GENOMES              } from '../../../modules/local/prepare_data_per_family'

include { paramsSummaryMap                                        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                    } from '../../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                  } from '../../nf-core/utils_nfcore_pipeline'
include { formatVersionsToYAML                                    } from '../utils_nfcore_chimeradetector_pipeline'
include { methodsDescriptionText                                  } from '../utils_nfcore_chimeradetector_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def getBlastDb( process_name ) {
    return process_name.tokenize(':')[1].tokenize('_')[2]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow MULTIQC_WORKFLOW {

    take:
    ch_chimeras_csv
    ch_versions

    main:

    ch_multiqc_files = Channel.empty()

    // ------------------------------------------------------------------------------------
    // VERSIONS
    // ------------------------------------------------------------------------------------

    // Collate and save software versions obtained from topic channels
    // TODO: use the nf-core functions when they are adapted to channel topics

    // Collate and save software versions
    formatVersionsToYAML ( Channel.topic('versions') )
        .mix ( softwareVersionsToYAML( ch_versions ) ) // mix with versions obtained from emit outputs
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }


    // ------------------------------------------------------------------------------------
    // CONFIG
    // ------------------------------------------------------------------------------------

    ch_multiqc_config = Channel.fromPath( "$projectDir/assets/multiqc_config.yml", checkIfExists: true )

    summary_params = paramsSummaryMap( workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value( paramsSummaryMultiqc(summary_params) )

    ch_multiqc_custom_config = params.multiqc_config ?
                                    Channel.fromPath(params.multiqc_config, checkIfExists: true) :
                                    Channel.empty()

    ch_multiqc_logo = params.multiqc_logo ?
                        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
                        Channel.empty()

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
                                                file(params.multiqc_methods_description, checkIfExists: true) :
                                                file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    ch_methods_description = Channel.value( methodsDescriptionText(ch_multiqc_custom_methods_description) )

    // Adding metadata to MultiQC
    ch_multiqc_files = ch_multiqc_files
                            .mix( ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml') )
                            .mix( ch_collated_versions )
                            .mix( ch_methods_description.collectFile( name: 'methods_description_mqc.yaml', sort: true ) )

    // ------------------------------------------------------------------------------------
    // DISTRIBUTION OF SIZE OF DOWNLAODED GENOMES PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('downloaded_genome_size')
        .collectFile(
            name: 'downloaded_genome_size.csv',
            seed: "family,data",
            newLine: true,
            storeDir: "${params.outdir}/assemblies/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_downloaded_genome_size_file }

    PREPARE_DOWNLOADED_GENOME_SIZE ( ch_downloaded_genome_size_file )

    // ------------------------------------------------------------------------------------
    // DISTRIBUTION OF SIZE OF ASSEMBLED GENOMES PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('assembled_genome_size')
        .collectFile(
            name: 'assembled_genome_size.csv',
            seed: "family,data",
            newLine: true,
            storeDir: "${params.outdir}/assemblies/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_assembled_genome_size_file }

    PREPARE_ASSEMBLED_GENOME_SIZE ( ch_assembled_genome_size_file )

    // ------------------------------------------------------------------------------------
    // DISTRIBUTION OF NB OF BLAST HITS AGAINST TARGET PER FAMILY
    // ------------------------------------------------------------------------------------

    Channel.topic('blast_nb_hits')
        .branch {
            process, family, nb_hits ->
                target: getBlastDb(process) == 'TARGET'
                    [ family, nb_hits ]
                genomes: getBlastDb(process) == 'GENOMES'
                   [ family, nb_hits ]
        }
        .set { ch_blast_nb_hits }

    ch_blast_nb_hits.target
        .collectFile(
            name: 'blast_nb_hits_target.csv',
            seed: "family,data",
            newLine: true,
            storeDir: "${params.outdir}/blastn/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_blast_nb_hits_target_file }

    PREPARE_BLAST_HIT_TARGET ( ch_blast_nb_hits_target_file )

    // ------------------------------------------------------------------------------------
    // DISTRIBUTION OF NB OF BLAST HITS AGAINST GENOMES PER FAMILY
    // ------------------------------------------------------------------------------------

    ch_blast_nb_hits.genomes
        .collectFile(
            name: 'blast_nb_hits_genomes.csv',
            seed: "family,data",
            newLine: true,
            storeDir: "${params.outdir}/blastn/"
        ) {
            item -> "${item[0]},${item[1]}"
        }
        .set { ch_blast_nb_hits_genomes_file }

    PREPARE_BLAST_HIT_GENOMES ( ch_blast_nb_hits_genomes_file )


    // ------------------------------------------------------------------------------------
    // LAUNCH MULTIQC
    // ------------------------------------------------------------------------------------

    ch_multiqc_files
        .mix ( ch_chimeras_csv )
        .mix ( PREPARE_DOWNLOADED_GENOME_SIZE.out.csv )
        .mix ( PREPARE_ASSEMBLED_GENOME_SIZE.out.csv )
        .mix ( PREPARE_BLAST_HIT_TARGET.out.csv )
        .mix ( PREPARE_BLAST_HIT_GENOMES.out.csv )
        .set { ch_multiqc_files }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )


    emit:
    multiqc_report = MULTIQC.out.report
}

