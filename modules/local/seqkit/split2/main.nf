process SEQKIT_SPLIT2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0' :
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads), val(nb_chunks)

    output:
    tuple val(meta), path("**/*.gz"), emit: reads
    tuple val("${task.process}"), val('seqkit'), eval("seqkit | sed '3!d; s/Version: //'"),     topic: versions

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqkit \\
        split2 \\
        $args \\
        --by-part $nb_chunks \\
        --threads $task.cpus \\
        $reads \\
        --out-dir ${prefix}
    """
}
