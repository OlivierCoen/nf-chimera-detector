process FIND_CHIMERAS {

    label 'process_low'

    tag "${meta.family} :: txid${meta.taxid} :: ${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7e/7e34603d81980331f4fb0224349f0af3019e234dc8c43130f9e3e36fff12dcd2/data':
        'community.wave.seqera.io/library/r-arrow_r-base_r-dplyr_r-optparse:55d208ef780e672c' }"

    input:
    tuple val(meta), path("blast_hits.against_target.txt"), path("blast_hits.against_genome.txt")

    output:
    tuple val(meta), path("*_chimeras.csv"),                                                                                  emit: csv
    tuple val("${meta.family}"),  val("${meta.id}"), env("NB_CHIMERAS"),                                                      topic: nb_chimeras
    tuple val("${task.process}"), val('R'),          eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'),       topic: versions
    tuple val("${task.process}"), val('dplyr'),      eval('Rscript -e "cat(as.character(packageVersion(\'dplyr\')))"'),       topic: versions
    tuple val("${task.process}"), val('data.table'), eval('Rscript -e "cat(as.character(packageVersion(\'data.table\')))"'),  topic: versions
    tuple val("${task.process}"), val('optparse'),   eval('Rscript -e "cat(as.character(packageVersion(\'optparse\')))"'),    topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_chimeras"
    """
    find_chimeras.R \\
        --hits-1 blast_hits.against_target.txt \\
        --hits-2 blast_hits.against_genome.txt \\
        --family ${meta.family} \\
        --species ${meta.taxid} \\
        --srr ${meta.id} \\
        --out ${prefix}.csv

    NB_CHIMERAS=\$(wc -l < chimeric_reads.txt)
    """

}
