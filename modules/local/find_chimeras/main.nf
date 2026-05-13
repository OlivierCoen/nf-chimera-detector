process FIND_CHIMERAS {

    label 'process_high'

    tag "${meta.family} :: txid${meta.taxid} :: ${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c4/c4f9fed0b61f39af8198f97aafc7c5cae07ab08570a44ebc48985995357628ce/data':
        'community.wave.seqera.io/library/r-arrow_r-base_r-data.table_r-dplyr_r-optparse:bcfb3d73ca929439' }"

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
